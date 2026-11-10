classdef DicomTree
    %DICOMTREE  Inventory of a heterogeneous DICOM tree: one row per series, classified by kind.
    %   t = mlphibroker.DicomTree(root) uses dicomCollection (Image Processing Toolbox) over all
    %   subfolders, then reads ONE header per series (dicominfo, header only) for PatientID,
    %   FrameOfReferenceUID, SeriesNumber and the PET frame count. Files that dicomCollection does
    %   not group (RTSTRUCT and other non-image objects) are swept separately with isdicom.
    %   Kinds: pet_dynamic, pet_static, ct, rtstruct, seg, other. PHI (names, dates, ids) stays
    %   in the in-memory table and is never written by this class.

    properties (SetAccess = private)
        Root
        Series table   % PatientID, StudyInstanceUID, SeriesInstanceUID, FrameOfReferenceUID,
                       % Modality, Kind, SeriesNumber, Frames, Files (cellstr), FirstHeader (struct)
    end

    methods
        function this = DicomTree(root)
            arguments
                root {mustBeFolder}
            end
            this.Root = char(root);
            coll = dicomCollection(this.Root, 'IncludeSubfolders', true);
            rows = {};
            seen = {};
            for i = 1:height(coll)
                files = cellstr(coll.Filenames{i});
                seen = [seen; files(:)]; %#ok<AGROW>
                rows{end+1} = mlphibroker.DicomTree.describe(files); %#ok<AGROW>
            end
            % Sweep for DICOM objects outside the collection (RTSTRUCT, SEG, SR, REG ...)
            all = dir(fullfile(this.Root, '**', '*'));
            all = all(~[all.isdir]);
            rest = setdiff(fullfile({all.folder}, {all.name}), seen);
            bySeries = containers.Map('KeyType', 'char', 'ValueType', 'any');
            for i = 1:numel(rest)
                p = rest{i};
                try
                    if ~isdicom(p), continue; end
                    h = dicominfo(p, 'UseDictionaryVR', true);
                catch
                    continue
                end
                if ~isfield(h, 'SeriesInstanceUID'), continue; end
                k = h.SeriesInstanceUID;
                if isKey(bySeries, k), v = bySeries(k); v{end+1} = p; bySeries(k) = v; %#ok<AGROW>
                else, bySeries(k) = {p}; end
            end
            ks = keys(bySeries);
            for i = 1:numel(ks)
                rows{end+1} = mlphibroker.DicomTree.describe(bySeries(ks{i})); %#ok<AGROW>
            end
            if isempty(rows)
                this.Series = table();
            else
                this.Series = vertcat(rows{:});
            end
        end

        function rows = seriesOfStudy(this, studyUid)
            rows = this.Series(strcmp(this.Series.StudyInstanceUID, studyUid), :);
        end

        function uids = orderStudiesByAcquisition(this, uids)
            % Sort study UIDs by their earliest StudyDate+StudyTime (else AcquisitionTime) read from
            % the first header of each series. PHI-side; only the resulting order is used.
            uids = cellstr(uids);
            if numel(uids) < 2, return; end
            key = zeros(numel(uids), 1);
            for i = 1:numel(uids)
                rows = this.seriesOfStudy(uids{i});
                best = inf;
                for r = 1:height(rows)
                    h = rows.FirstHeader{r};
                    d = ''; t = '';
                    if isfield(h, 'StudyDate'), d = char(string(h.StudyDate)); end
                    if isfield(h, 'SeriesTime'), t = char(string(h.SeriesTime)); elseif isfield(h, 'AcquisitionTime'), t = char(string(h.AcquisitionTime)); end
                    v = str2double([regexprep(d, '\D', '') sprintf('%06.0f', mlphibroker.PetSession.clockSeconds(t))]);
                    if ~isnan(v), best = min(best, v); end
                end
                key(i) = best;
            end
            [~, order] = sort(key);
            uids = uids(order);
        end

        function c = countByKind(this)
            c = struct();
            kinds = {'pet_dynamic', 'pet_static', 'ct', 'rtstruct', 'seg', 'secondary', 'localizer', 'other'};
            for i = 1:numel(kinds)
                c.(kinds{i}) = sum(strcmp(this.Series.Kind, kinds{i}));
            end
        end

        function runOf = runNumbers(this, studyUids)
            %RUNNUMBERS  Map SeriesInstanceUID -> run index per kind across the given studies (0 = unnumbered).
            %   Studies are taken in the order given; within a study, series are ordered by their first
            %   header's SeriesTime/AcquisitionTime (PHI-side, only the order is kept), then SeriesNumber.
            %   A kind with a single series in the session stays unnumbered (plain file names).
            runOf = containers.Map('KeyType', 'char', 'ValueType', 'double');
            for kindSet = {{'pet_dynamic', 'pet_static'}, {'ct'}}
                uids = {}; keys = [];
                for st = 1:numel(studyUids)
                    rows = this.seriesOfStudy(studyUids{st});
                    rows = rows(ismember(rows.Kind, kindSet{1}), :);
                    for r = 1:height(rows)
                        h = rows.FirstHeader{r}; t = '';
                        if isfield(h, 'SeriesTime'), t = char(string(h.SeriesTime)); elseif isfield(h, 'AcquisitionTime'), t = char(string(h.AcquisitionTime)); end
                        sec = mlphibroker.PetSession.clockSeconds(t); if isempty(sec), sec = 0; end
                        uids{end+1} = rows.SeriesInstanceUID{r}; %#ok<AGROW>
                        keys(end+1) = st * 1e6 + sec + rows.SeriesNumber(r) * 1e-6; %#ok<AGROW>
                    end
                end
                [~, order] = sort(keys);
                for i = 1:numel(order)
                    if numel(order) > 1, runOf(uids{order(i)}) = i; else, runOf(uids{order(i)}) = 0; end
                end
            end
            for st = 1:numel(studyUids)   % every series gets an entry so lookups never miss
                rows = this.seriesOfStudy(studyUids{st});
                for r = 1:height(rows)
                    if ~isKey(runOf, rows.SeriesInstanceUID{r}), runOf(rows.SeriesInstanceUID{r}) = 0; end
                end
            end
        end

        function c = countUnhandledByModality(this)
            c = struct();
            rows = this.Series(ismember(this.Series.Kind, {'other', 'secondary', 'localizer'}), :);
            for m = unique(rows.Modality)'
                key = regexprep(m{1}, '[^A-Za-z0-9]', '_'); if isempty(key), key = 'blank'; end
                c.(key) = sum(strcmp(rows.Modality, m{1}));
            end
        end

        function ids = uniquePatients(this)
            ids = unique(this.Series.PatientID, 'stable');
        end
    end

    methods (Static)
        function inspectNames(root, kinds)
            %INSPECTNAMES  PHI-SIDE CONSOLE ONLY: print every name-bearing or free-text header field of
            %   one series per kind, so the operator can decide what the deny list must contain.
            %   mlphibroker.DicomTree.inspectNames(root)            PET (dynamic, static), CT, RTSTRUCT
            %   mlphibroker.DicomTree.inspectNames(root, {'rtstruct'})
            %   Prints: every field whose value is a DICOM Person Name (any tag), the usual free-text
            %   fields, and for an RTSTRUCT the structure-set labels, ROI names and ROI interpreters.
            %   The output shows values: never paste it outside the PHI zone.
            if nargin < 2 || isempty(kinds), kinds = {'pet_dynamic', 'pet_static', 'ct', 'rtstruct'}; end
            tree = mlphibroker.DicomTree(root);
            free = {'StudyDescription', 'SeriesDescription', 'ProtocolName', 'ImageComments', ...
                'RequestedProcedureDescription', 'PerformedProcedureStepDescription', 'InstitutionName', ...
                'InstitutionalDepartmentName', 'InstitutionAddress', 'StationName', 'DeviceSerialNumber', ...
                'PatientID', 'OtherPatientIDs', 'AccessionNumber', 'StudyID', 'PatientComments', ...
                'AdditionalPatientHistory', 'PatientBirthDate', 'PatientAge', 'PatientSex', ...
                'StructureSetLabel', 'StructureSetName', 'StructureSetDescription'};
            for kind = cellstr(kinds)
                rows = tree.Series(strcmp(tree.Series.Kind, kind{1}), :);
                if isempty(rows), fprintf('\n== %s: none in tree\n', kind{1}); continue; end
                for r = 1:height(rows)
                    h = rows.FirstHeader{r};
                    fprintf('\n== %s  ser-%g  (%d files)\n', kind{1}, rows.SeriesNumber(r), numel(rows.Files{r}));
                    f = fieldnames(h);
                    for i = 1:numel(f)
                        v = h.(f{i});
                        if isstruct(v) && isfield(v, 'FamilyName')
                            fprintf('  PN  %-34s %s\n', f{i}, mlphibroker.DicomTree.pnText(v));
                        end
                    end
                    for k = free
                        if isfield(h, k{1}) && ~isempty(h.(k{1}))
                            fprintf('  TXT %-34s %s\n', k{1}, mlphibroker.DicomTree.anyText(h.(k{1})));
                        end
                    end
                    if strcmp(kind{1}, 'rtstruct')
                        if isfield(h, 'StructureSetROISequence')
                            q = h.StructureSetROISequence;
                            for n = fieldnames(q)'
                                it = q.(n{1});
                                fprintf('  ROI %-34s %s\n', 'ROIName', mlphibroker.DicomTree.anyText(mlphibroker.DicomTree.field(it, 'ROIName', '')));
                            end
                        end
                        if isfield(h, 'RTROIObservationsSequence')
                            q = h.RTROIObservationsSequence;
                            for n = fieldnames(q)'
                                it = q.(n{1});
                                for k = {'ROIInterpreter', 'ROIObservationLabel', 'ROIObservationDescription'}
                                    v = mlphibroker.DicomTree.field(it, k{1}, '');
                                    if ~isempty(v), fprintf('  ROI %-34s %s\n', k{1}, mlphibroker.DicomTree.anyText(v)); end
                                end
                            end
                        end
                    end
                end
            end
        end

        function t = pnText(v)
            % Person Name struct -> 'Family^Given^Middle^Prefix^Suffix' with empty parts dropped.
            parts = {};
            for k = {'FamilyName', 'GivenName', 'MiddleName', 'NamePrefix', 'NameSuffix'}
                if isfield(v, k{1}) && ~isempty(v.(k{1})), parts{end+1} = char(string(v.(k{1}))); end %#ok<AGROW>
            end
            if isempty(parts), t = '(empty)'; else, t = strjoin(parts, '^'); end
        end

        function t = anyText(v)
            % Any header value -> one printable line (PN structs, numbers, arrays, cellstr).
            if isstruct(v) && isfield(v, 'FamilyName'), t = mlphibroker.DicomTree.pnText(v);
            elseif isstruct(v), t = sprintf('<struct with %d fields>', numel(fieldnames(v)));
            elseif iscell(v), t = strjoin(cellfun(@(x) mlphibroker.DicomTree.anyText(x), v, 'UniformOutput', false), ' | ');
            elseif isnumeric(v) || islogical(v), t = mat2str(double(v(:))', 6);
            else, t = char(string(v)); end
        end
    end

    methods (Static, Access = private)
        function row = describe(files)
            files = cellstr(files);
            h = dicominfo(files{1}, 'UseDictionaryVR', true);
            get = @(f, d) mlphibroker.DicomTree.field(h, f, d);
            modality = upper(char(get('Modality', '')));
            frames = 1;
            sop = char(get('SOPClassUID', ''));
            itype = get('ImageType', '');
            if ~ischar(itype), itype = strjoin(cellfun(@char, cellstr(string(itype)), 'UniformOutput', false), '\\'); end
            itype = upper(itype);
            isSecondary = startsWith(sop, '1.2.840.10008.5.1.4.1.1.7') || contains(itype, 'SECONDARY') || ...
                contains(itype, 'SCREEN') || contains(itype, 'DOSE_INFO') || contains(itype, 'MIP');
            isLocalizer = contains(itype, 'LOCALIZER');
            if isSecondary, modality = ['SECONDARY(' modality ')']; end
            if isLocalizer, modality = ['LOCALIZER(' modality ')']; end
            switch modality
                case 'PT'
                    nts = get('NumberOfTimeSlices', []);
                    if ~isempty(nts) && nts > 1
                        frames = double(nts);
                    else
                        frames = mlphibroker.DicomTree.countAcquisitions(files);
                    end
                    if frames > 1, kind = 'pet_dynamic'; else, kind = 'pet_static'; end
                case 'CT', kind = 'ct';
                case 'RTSTRUCT', kind = 'rtstruct';
                case 'SEG', kind = 'seg';
                otherwise
                    if isSecondary, kind = 'secondary';       % screen saves, dose reports, MIPs: text may be burned into pixels
                    elseif isLocalizer, kind = 'localizer';   % topograms: annotations possible, no quantitative use
                    else, kind = 'other'; end
            end
            row = table({char(get('PatientID', ''))}, {char(get('StudyInstanceUID', ''))}, ...
                {char(get('SeriesInstanceUID', ''))}, {char(get('FrameOfReferenceUID', ''))}, ...
                {modality}, {kind}, double(get('SeriesNumber', NaN)), frames, {files(:)'}, {h}, ...
                'VariableNames', {'PatientID', 'StudyInstanceUID', 'SeriesInstanceUID', ...
                'FrameOfReferenceUID', 'Modality', 'Kind', 'SeriesNumber', 'Frames', 'Files', 'FirstHeader'});
        end

        function v = field(h, f, d)
            if isfield(h, f) && ~isempty(h.(f)), v = h.(f); else, v = d; end
        end

        function n = countAcquisitions(files)
            % Distinct AcquisitionNumber over at most 64 sampled headers (cheap dynamic test).
            idx = unique(round(linspace(1, numel(files), min(64, numel(files)))));
            acq = [];
            for i = idx
                try
                    h = dicominfo(files{i}, 'UseDictionaryVR', true);
                    if isfield(h, 'AcquisitionNumber'), acq(end+1) = double(h.AcquisitionNumber); end %#ok<AGROW>
                catch
                end
            end
            n = max(1, numel(unique(acq)));
        end
    end
end
