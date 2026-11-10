classdef Crosswalk
    %CROSSWALK  PatientID -> sub code and StudyInstanceUID -> random ses ordinal (PHI-side).
    %   xw = mlphibroker.Crosswalk.resolve(tree, outRoot, opts) returns a table with columns
    %   PatientID, Subject, StudyInstanceUID, Session, one row per study in the tree.
    %   Source: opts.Crosswalk CSV from the key vault; else, if opts.AutoCrosswalk, synthetic
    %   codes are minted (sub-9000NN in first-appearance order; ses-NN a random permutation per
    %   subject, SPEC-PHIBROKER §4.2) and written to outRoot/crosswalk_PHI.csv, which holds
    %   PatientIDs and MUST stay in the PHI zone (.gitignore excludes *_PHI.csv).
    %   A study in the tree with no crosswalk row is an error: codes are never guessed.

    methods (Static)
        function xw = resolve(tree, outRoot, opts)
            studies = unique(tree.Series(:, {'PatientID', 'StudyInstanceUID'}), 'stable');
            if isfield(opts, 'Crosswalk') && ~isempty(opts.Crosswalk)
                xw = mlphibroker.Crosswalk.load(opts.Crosswalk);
            elseif isfield(opts, 'AutoCrosswalk') && isequal(opts.AutoCrosswalk, true)
                xw = mlphibroker.Crosswalk.mint(studies, opts);
                writetable(xw, fullfile(outRoot, 'crosswalk_PHI.csv'));
            else
                error('mlphibroker:crosswalk', ...
                    'opts.Crosswalk (key-vault CSV) is required, or opts.AutoCrosswalk=true for synthetic codes');
            end
            mlphibroker.Crosswalk.validate(xw, mlphibroker.Crosswalk.subjectPattern(opts));
            missing = ~ismember(studies.StudyInstanceUID, xw.StudyInstanceUID);
            if any(missing)
                error('mlphibroker:crosswalk', '%d studies in the tree have no crosswalk row', sum(missing));
            end
            xw = xw(ismember(xw.StudyInstanceUID, studies.StudyInstanceUID), :);
        end

        function t = template(root, outPath, opts)
            %TEMPLATE  Inventory a DICOM tree and write a PHI-side crosswalk with PROPOSED codes to authorize.
            %   t = mlphibroker.Crosswalk.template(root, outPath, opts) writes one row per study:
            %     PatientID, PatientName, StudyDate, StudyDescription, Modalities, StudyInstanceUID
            %       (read from the DICOM headers; PHI, for your eyes only)
            %     Subject   proposed: opts.SubjectPrefix + 3-digit number in first-appearance order, one
            %               per PERSON, where a person is one PatientName even under several PatientIDs
            %               (MIM exports often split them). Empty if opts.SubjectPrefix is not given.
            %     Session   proposed: ses-NN. Studies of one person on one StudyDate are ONE visit and
            %               share a Session; distinct days get a RANDOM ordinal (SPEC-PHIBROKER §4.2), so
            %               ses-01 is not necessarily the earliest visit. Seed: opts.Seed (default 0).
            %     Note      why each value was proposed. Edit Subject/Session if you disagree; everything
            %               else is read-only context. load() keeps only the four crosswalk columns.
            %   A companion <outPath>.README.txt repeats this. The CSV holds PHI: keep it in the PHI zone
            %   and name it *_PHI.csv so the repository's .gitignore excludes it.
            if nargin < 3 || ~isstruct(opts), opts = struct(); end
            tree = mlphibroker.DicomTree(root);
            S = tree.Series;
            [~, first] = unique(S.StudyInstanceUID, 'stable');
            n = numel(first);
            PatientID = cell(n, 1); PatientName = cell(n, 1); StudyDate = cell(n, 1);
            StudyDescription = cell(n, 1); StudyInstanceUID = cell(n, 1); Modalities = cell(n, 1);
            for i = 1:n
                h = S.FirstHeader{first(i)};
                PatientID{i} = S.PatientID{first(i)};
                PatientName{i} = mlphibroker.Crosswalk.nameOf(h);
                StudyDate{i} = mlphibroker.Crosswalk.strField(h, 'StudyDate');
                StudyDescription{i} = mlphibroker.Crosswalk.strField(h, 'StudyDescription');
                StudyInstanceUID{i} = S.StudyInstanceUID{first(i)};
                rows = S(strcmp(S.StudyInstanceUID, StudyInstanceUID{i}), :);
                Modalities{i} = strjoin(unique(rows.Modality)', '+');
            end
            % Person key: the name when present (ids may differ per export), else the id.
            person = cellfun(@(nm, id) lower(strtrim(iff(~isempty(nm), nm, id))), PatientName, PatientID, 'UniformOutput', false);
            [persons, ~, pidx] = unique(person, 'stable');
            Subject = repmat({''}, n, 1); Session = repmat({''}, n, 1); Note = repmat({''}, n, 1);
            prefix = ''; if isfield(opts, 'SubjectPrefix'), prefix = char(opts.SubjectPrefix); end
            seed = 0; if isfield(opts, 'Seed'), seed = opts.Seed; end
            rs = RandStream('mt19937ar', 'Seed', seed);
            roster = table(); if isfield(opts, 'Roster') && ~isempty(opts.Roster), roster = mlphibroker.Crosswalk.loadRoster(opts.Roster); end
            for k = 1:numel(persons)
                rowsK = find(pidx == k);
                ids = unique(PatientID(rowsK));
                days = unique(StudyDate(rowsK), 'stable');
                perm = randperm(rs, numel(days));
                for r = rowsK'
                    note = {};
                    code = '';
                    if ~isempty(roster)
                        code = mlphibroker.Crosswalk.rosterCode(roster, PatientID{r}, PatientName{r});
                    end
                    if ~isempty(code)
                        Subject{r} = code;
                        note{end+1} = 'Subject from the coordinator roster (authoritative)'; %#ok<AGROW>
                    elseif ~isempty(roster)
                        note{end+1} = 'NOT IN ROSTER: ask the study coordinator for this person''s code; left empty'; %#ok<AGROW>
                    elseif ~isempty(prefix)
                        Subject{r} = sprintf('%s%03d', prefix, k);
                        note{end+1} = sprintf('PLACEHOLDER %d of %d by first appearance in this export; replace with the coordinator''s code', k, numel(persons)); %#ok<AGROW>
                    else
                        note{end+1} = 'Subject: give opts.Roster (coordinator codes) or opts.SubjectPrefix (placeholders)'; %#ok<AGROW>
                    end
                    if numel(ids) > 1
                        note{end+1} = sprintf('%d PatientIDs share this name: same Subject proposed', numel(ids)); %#ok<AGROW>
                    end
                    d = find(strcmp(days, StudyDate{r}), 1);
                    Session{r} = sprintf('ses-%02d', perm(d));
                    sameDay = sum(strcmp(StudyDate(rowsK), StudyDate{r}));
                    if sameDay > 1
                        note{end+1} = sprintf('%d studies on one day = one visit: same Session (outputs run-01, run-02 by acquisition order)', sameDay); %#ok<AGROW>
                    end
                    if numel(days) > 1
                        note{end+1} = sprintf('%d visit days: Session ordinals assigned at random, not by date', numel(days)); %#ok<AGROW>
                    end
                    Note{r} = strjoin(note, '; ');
                end
            end
            t = table(PatientID, PatientName, StudyDate, StudyDescription, Modalities, ...
                StudyInstanceUID, Subject, Session, Note);
            t = sortrows(t, {'PatientName', 'StudyDate'});
            writetable(t, outPath);
            mlphibroker.Crosswalk.writeTemplateReadme([outPath '.README.txt'], prefix, seed);
            fprintf('crosswalk template: %d studies, %d persons (%d PatientIDs) -> %s (PHI; authorize or edit Subject and Session; see %s)\n', ...
                n, numel(persons), numel(unique(PatientID)), outPath, [outPath '.README.txt']);

            function v = iff(c, a, b)
                if c, v = a; else, v = b; end
            end
        end

        function roster = loadRoster(path)
            %LOADROSTER  Coordinator roster: columns Subject and PatientID and/or PatientName (PHI-side).
            if ~isfile(path), error('mlphibroker:roster', 'roster file not found'); end
            roster = readtable(path, 'TextType', 'char', 'Delimiter', ',');
            if ~ismember('Subject', roster.Properties.VariableNames) || ...
                    ~any(ismember({'PatientID', 'PatientName'}, roster.Properties.VariableNames))
                error('mlphibroker:roster', 'roster needs a Subject column and a PatientID and/or PatientName column');
            end
            for f = intersect({'Subject', 'PatientID', 'PatientName'}, roster.Properties.VariableNames)
                roster.(f{1}) = cellstr(string(roster.(f{1})));
            end
        end

        function code = rosterCode(roster, pid, pname)
            % Match by PatientID first, then by case-insensitive PatientName (DICOM caret form or spaces).
            code = '';
            hit = false(height(roster), 1);
            if ismember('PatientID', roster.Properties.VariableNames) && ~isempty(pid)
                hit = strcmp(strtrim(roster.PatientID), strtrim(pid));
            end
            if ~any(hit) && ismember('PatientName', roster.Properties.VariableNames) && ~isempty(pname)
                norm = @(x) lower(regexprep(strtrim(x), '[\s\^,]+', '^'));
                hit = strcmp(cellfun(norm, roster.PatientName, 'UniformOutput', false), norm(pname));
            end
            codes = unique(roster.Subject(hit));
            if numel(codes) == 1, code = codes{1}; end
        end

        function writeTemplateReadme(path, prefix, seed)
            lines = {
                'mlphibroker crosswalk template: how to complete it'
                ''
                'This CSV holds PHI (names, ids, dates). Keep it in the PHI zone; never commit or share it.'
                'Only the columns PatientID, Subject, StudyInstanceUID, Session are used by the conversion.'
                ''
                'Columns read from the DICOM headers (context, do not edit):'
                '  PatientID, PatientName, StudyDate, StudyDescription, Modalities, StudyInstanceUID'
                ''
                'Columns to AUTHORIZE (proposed values are filled in; edit only if you disagree):'
                '  Subject   the study coordinator''s code for the person. With opts.Roster (a CSV of Subject and'
                '            PatientID and/or PatientName kept in the PHI zone) it is filled from the roster and is'
                '            authoritative; a person missing from the roster is left EMPTY and must be asked for.'
                sprintf('            Without a roster, opts.SubjectPrefix ("%s") fills PLACEHOLDERS numbered by first', prefix)
                '            appearance in this export; such numbers do NOT match the coordinator''s (contacts that'
                '            never enrolled leave gaps), so replace them before converting anything you will keep.'
                '            A person is one PatientName, even when the export used several PatientIDs.'
                '            Allowed shapes: sub-108xxx (CCIR_01211), sub-ECL1i001-style (MASH), sub-9xxxxx (synthetic);'
                '            labels are alphanumeric (sub-ECL1i-001 is refused).'
                '  Session   ses-NN per visit. Studies of one person on one StudyDate are ONE visit and share'
                '            a Session (their outputs get run-01, run-02 ... by acquisition order).'
                sprintf('            Distinct days get RANDOM ordinals (seed %d), so ses-01 is not necessarily the earliest visit;', seed)
                '            this is deliberate (SPEC-PHIBROKER section 4.2): the ordinal must not encode chronology.'
                '  Note      explains each proposal; not read by the conversion.'
                ''
                'Refused at load time: an empty Subject or Session; a hyphen inside a subject label; a Subject'
                'outside the allowed shapes (opts.SubjectPattern admits another project); a StudyInstanceUID'
                'listed twice; a study present in the DICOM tree but absent from this file.'
                };
            fid = fopen(path, 'w'); fprintf(fid, '%s\n', lines{:}); fclose(fid);
        end

        function xw = load(path)
            if ~isfile(path), error('mlphibroker:crosswalk', 'crosswalk file not found'); end
            xw = readtable(path, 'TextType', 'char', 'Delimiter', ',');
            need = {'PatientID', 'Subject', 'StudyInstanceUID', 'Session'};
            for i = 1:numel(need)
                if ~ismember(need{i}, xw.Properties.VariableNames)
                    error('mlphibroker:crosswalk', 'crosswalk lacks column %s', need{i});
                end
                xw.(need{i}) = cellstr(string(xw.(need{i})));
            end
            blank = cellfun(@(s) isempty(strtrim(s)) || strcmp(s, '<missing>'), xw.Subject) | ...
                cellfun(@(s) isempty(strtrim(s)) || strcmp(s, '<missing>'), xw.Session);
            if any(blank)
                error('mlphibroker:crosswalk', '%d crosswalk rows have an empty Subject or Session (fill the template first)', sum(blank));
            end
            xw = xw(:, need);   % extra template columns (names, dates) are never carried further
        end

        function xw = mint(studies, opts)
            pids = unique(studies.PatientID, 'stable');
            seed = 0; if isfield(opts, 'Seed'), seed = opts.Seed; end
            rs = RandStream('mt19937ar', 'Seed', seed);
            Subject = cell(height(studies), 1); Session = cell(height(studies), 1);
            for i = 1:numel(pids)
                idx = find(strcmp(studies.PatientID, pids{i}));
                perm = randperm(rs, numel(idx));
                for k = 1:numel(idx)
                    Subject{idx(k)} = sprintf('sub-9%05d', i);
                    Session{idx(k)} = sprintf('ses-%02d', perm(k));
                end
            end
            xw = table(studies.PatientID, Subject, studies.StudyInstanceUID, Session, ...
                'VariableNames', {'PatientID', 'Subject', 'StudyInstanceUID', 'Session'});
        end

        function v = strField(h, f)
            if isfield(h, f) && ~isempty(h.(f)), v = char(string(h.(f))); else, v = ''; end
        end

        function v = nameOf(h)
            v = '';
            if isfield(h, 'PatientName')
                pn = h.PatientName;
                if isstruct(pn)
                    parts = {};
                    for f = {'FamilyName', 'GivenName'}
                        if isfield(pn, f{1}) && ~isempty(pn.(f{1})), parts{end+1} = char(pn.(f{1})); end %#ok<AGROW>
                    end
                    v = strjoin(parts, '^');
                else
                    v = char(string(pn));
                end
            end
        end

        function pat = subjectPattern(opts)
            %SUBJECTPATTERN  Accepted subject codes: CCIR_01211 sub-108xxx, synthetic sub-9xxxxx, and
            %   the MASH (DOTA-ECL1i) range sub-ECL1i followed by three digits (operator, 2026-11-10).
            %   opts.SubjectPattern (a regex) overrides for another project. Labels stay BIDS
            %   alphanumerics: a hyphen inside the label (sub-ECL1i-001) breaks sub-/ses-/trc- parsing.
            pat = '^sub-(108\d{3}|9\d{5}|ECL1i\d{3})$';
            if nargin > 0 && isstruct(opts) && isfield(opts, 'SubjectPattern') && ~isempty(opts.SubjectPattern)
                pat = char(opts.SubjectPattern);
            end
        end

        function validate(xw, pat)
            if nargin < 2 || isempty(pat), pat = mlphibroker.Crosswalk.subjectPattern(struct()); end
            hyph = cellfun(@(s) ~isempty(regexp(s, '^sub-[A-Za-z0-9]+-', 'once')), xw.Subject);
            if any(hyph)
                error('mlphibroker:crosswalk', ['%d Subject codes carry a hyphen inside the label (e.g. sub-ECL1i-001); ' ...
                    'BIDS labels are alphanumeric: use sub-ECL1i001'], sum(hyph));
            end
            bad = cellfun(@(s) isempty(regexp(s, pat, 'once')), xw.Subject);
            if any(bad), error('mlphibroker:crosswalk', '%d Subject codes do not match %s', sum(bad), pat); end
            bad = cellfun(@(s) isempty(regexp(s, '^ses-\d{2}$', 'once')), xw.Session);
            if any(bad), error('mlphibroker:crosswalk', '%d Session codes are not ses-NN (dated labels are refused)', sum(bad)); end
            % Several rows may share Subject/Session: one visit with several studies (e.g. a dynamic
            % and a delayed scan of one injection); their outputs are told apart by run-NN. A
            % StudyInstanceUID may appear once only.
            if numel(unique(xw.StudyInstanceUID)) ~= numel(xw.StudyInstanceUID)
                error('mlphibroker:crosswalk', 'a StudyInstanceUID appears more than once in the crosswalk');
            end
            pid = cell(0, 1);
            for sub = unique(xw.Subject)'
                pid = xw.PatientID(strcmp(xw.Subject, sub{1}));
                if numel(unique(pid)) > 1
                    fprintf('note: %s maps %d PatientIDs to one subject (accepted: same person, several ids)\n', sub{1}, numel(unique(pid)));
                end
            end
        end
    end
end
