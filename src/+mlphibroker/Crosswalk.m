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

        function t = template(root, outPath)
            %TEMPLATE  Inventory a DICOM tree and write a PHI-side crosswalk template to fill in.
            %   t = mlphibroker.Crosswalk.template(root, outPath) writes one row per study with
            %   PatientID, PatientName, StudyDate, StudyDescription, StudyInstanceUID, Modalities,
            %   and EMPTY Subject and Session columns. The operator fills Subject (sub-108xxx, or
            %   sub-9xxxxx for synthetic data) and Session (ses-NN, a RANDOM ordinal per subject,
            %   never derived from the date), then passes the file as opts.Crosswalk.
            %   The template holds PHI (names, ids, dates): it stays in the PHI zone. Name it
            %   *_PHI.csv so the repository's .gitignore excludes it.
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
            Subject = repmat({''}, n, 1); Session = repmat({''}, n, 1);
            t = table(PatientID, PatientName, StudyDate, StudyDescription, Modalities, ...
                StudyInstanceUID, Subject, Session);
            t = sortrows(t, {'PatientName', 'StudyDate'});
            writetable(t, outPath);
            [~, ~, gi] = unique(t.PatientName);
            split = 0;
            for g = 1:max(gi), if numel(unique(t.PatientID(gi == g))) > 1, split = split + 1; end, end
            if split > 0
                fprintf(['note: %d patient name(s) carry more than one PatientID; if they are one person, ' ...
                    'give those rows the same Subject. Studies of one visit (e.g. dynamic + delayed) take the same Session.\n'], split);
            end
            fprintf('crosswalk template: %d studies, %d patients -> %s (PHI; fill Subject and Session)\n', ...
                n, numel(unique(PatientID)), outPath);
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
