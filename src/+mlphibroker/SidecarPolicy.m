classdef SidecarPolicy
    %SIDECARPOLICY  Allow-list for sidecar keys, deny-list of names, and the PHI-shape scan.
    %   p = mlphibroker.SidecarPolicy(denyNames)
    %   [kept, dropped] = p.allow(sidecarStruct)      keep listed keys only (fail-closed)
    %   hits = p.scan(text)                            counts of date / clock / dated-ses / deny hits
    %   Mirrors the $Allow table of idif2026 scripts/phibroker/enclave/Invoke-Dcm2niixDeid.ps1
    %   (WO-0002) plus dcm2niix spellings observed in real output: BodyPart, ReconstructionMethod
    %   (v1.0.20250505) and the BIDS-PET names TracerName, DecayCorrectionFactor (v1.0.20261029),
    %   and CT acquisition parameters. dcm2niix 2026 also writes SeriesTime, a clock time, which
    %   this list drops like AcquisitionTime and TimeZero. Identity, every absolute date or clock time
    %   (AcquisitionTime, TimeZero, ...), serials, station and institution names and UIDs are
    %   never on the list. Changing the list is a reviewed policy change (OQ-14).

    properties (Constant)
        ALLOW = { ...
            'Modality', 'Manufacturer', 'ManufacturersModelName', 'SoftwareVersions', ...
            'BodyPartExamined', 'BodyPart', 'PatientPosition', ...
            'SeriesDescription', 'ProtocolName', 'ImageType', 'SeriesNumber', 'AcquisitionNumber', ...
            'SliceThickness', 'PixelSpacing', 'ImageOrientationPatientDICOM', ...
            'Units', 'ReconMethodName', 'ReconstructionMethod', 'ReconMethodParameterLabels', ...
            'ReconMethodParameterUnits', 'ReconMethodParameterValues', 'ReconFilterType', ...
            'ReconFilterSize', 'AttenuationCorrection', 'AttenuationCorrectionMethod', ...
            'ScatterCorrectionMethod', 'RandomsCorrectionMethod', 'ScatterFraction', ...
            'DoseCalibrationFactor', ...
            'DecayCorrection', 'DecayFactor', 'DecayCorrectionFactor', 'ImageDecayCorrected', 'ImageDecayCorrectionTime', ...
            'Radiopharmaceutical', 'TracerName', 'TracerRadionuclide', 'RadionuclideTotalDose', 'RadionuclideHalfLife', ...
            'RadionuclidePositronFraction', 'InjectedRadioactivity', 'InjectedRadioactivityUnits', ...
            'SpecificRadioactivity', 'SpecificRadioactivityUnits', 'ModeOfAdministration', ...
            'FrameTimesStart', 'FrameDuration', 'FrameReferenceTime', 'ScanStart', 'InjectionStart', ...
            'PatientWeight', ...
            'KVP', 'XRayExposure', 'XRayTubeCurrent', 'ConvolutionKernel', 'ExposureTime', ...
            'ReconstructionDiameter', 'SpiralPitchFactor', 'CTDIvol', 'TotalCollimationWidth', ...
            'SingleCollimationWidth', 'TableSpeed', 'TableFeedPerRotation', ...
            'ConversionSoftware', 'ConversionSoftwareVersion'}
    end

    properties (SetAccess = private)
        Deny cell = {}
    end

    methods
        function this = SidecarPolicy(deny)
            if nargin > 0, this.Deny = cellstr(deny); end
        end

        function [kept, dropped] = allow(this, src)
            kept = struct(); dropped = {};
            f = fieldnames(src);
            for i = 1:numel(f)
                if any(strcmp(f{i}, this.ALLOW)), kept.(f{i}) = src.(f{i});
                else, dropped{end+1} = f{i}; end %#ok<AGROW>
            end
        end

        function hits = scan(this, txt)
            txt = char(txt);
            hits = struct('date', 0, 'clock_time', 0, 'dated_session', 0, 'deny_list', 0, 'any', false);
            hits.date = numel(regexp(txt, '(?<![0-9])(19|20)\d{2}-[01]\d-[0-3]\d(?![0-9])')) ...
                + numel(regexp(txt, '(?<![0-9.])(19|20)\d{2}[01]\d[0-3]\d(?![0-9.])')) ...
                + numel(regexp(txt, '(?<![0-9])[01]?\d/[0-3]?\d/(19|20)\d{2}(?![0-9])'));
            hits.clock_time = numel(regexp(txt, '(?<![0-9])[0-2]\d:[0-5]\d(:[0-5]\d(\.\d+)?)?(?![0-9])'));
            hits.dated_session = numel(regexp(txt, 'ses-\d{8,}'));
            for i = 1:numel(this.Deny)
                d = strtrim(char(this.Deny{i}));
                if numel(d) < 3, continue; end
                hits.deny_list = hits.deny_list + numel(regexpi(txt, ...
                    ['(?<![A-Za-z])' regexptranslate('escape', d) '(?![A-Za-z])']));
            end
            hits.any = (hits.date + hits.clock_time + hits.dated_session + hits.deny_list) > 0;
        end

        function [hits, byField] = scanStruct(this, sc, extraText)
            %SCANSTRUCT  Scan a sidecar field by field; byField lists the KEYS that carry a hit.
            %   Key names are not PHI and may be reported; values never are.
            if nargin < 3, extraText = ''; end
            hits = this.scan([jsonencode(sc) ' ' extraText]);
            byField = {};
            f = fieldnames(sc);
            for i = 1:numel(f)
                h = this.scan(jsonencode(sc.(f{i})));
                if h.any
                    kinds = {};
                    if h.date > 0, kinds{end+1} = 'date'; end %#ok<AGROW>
                    if h.clock_time > 0, kinds{end+1} = 'clock_time'; end %#ok<AGROW>
                    if h.dated_session > 0, kinds{end+1} = 'dated_session'; end %#ok<AGROW>
                    if h.deny_list > 0, kinds{end+1} = 'deny_list'; end %#ok<AGROW>
                    byField{end+1} = sprintf('%s:%s', f{i}, strjoin(kinds, '+')); %#ok<AGROW>
                end
            end
            if ~isempty(extraText) && this.scan(extraText).any, byField{end+1} = 'file_name'; end
        end
    end

    methods (Static)
        function deny = readDenyList(path)
            if ~isfile(path), error('mlphibroker:deny', 'DenyListFile not found'); end
            lines = regexp(fileread(path), '\r?\n', 'split');
            deny = {};
            for i = 1:numel(lines)
                d = strtrim(lines{i});
                if numel(d) >= 3, deny{end+1} = d; end %#ok<AGROW>
            end
        end

        function names = denyListTemplate(root, outPath, extra)
            %DENYLISTTEMPLATE  Build a deny list from every name-bearing header field in a DICOM tree.
            %   names = mlphibroker.SidecarPolicy.denyListTemplate(root, outPath, extra) reads one
            %   header per series (PatientName, OperatorsName, ReferringPhysicianName,
            %   PerformingPhysicianName, PhysiciansOfRecord, RequestingPhysician,
            %   NameOfPhysiciansReadingStudy, InstitutionName, InstitutionalDepartmentName,
            %   StationName, DeviceSerialNumber, PatientID, AccessionNumber, StudyID), splits names into words of three
            %   or more letters, keeps the full strings too, adds `extra` (cellstr: coordinator,
            %   investigators, anyone who might appear in free text) and writes one entry per line.
            %   PHI: the file stays in the PHI zone. Review it and add names before converting.
            if nargin < 3, extra = {}; end
            % Person-name fields are split into words (surname and given name match on their own);
            % identifiers, institution and station are kept as WHOLE strings only: splitting an id such
            % as ECL1i-005 would yield 'ECL' and quarantine every protocol named after the tracer.
            nameFields = {'PatientName', 'OperatorsName', 'ReferringPhysicianName', 'PerformingPhysicianName', ...
                'PhysiciansOfRecord', 'RequestingPhysician', 'NameOfPhysiciansReadingStudy'};
            wholeFields = {'InstitutionName', 'InstitutionalDepartmentName', 'StationName', 'DeviceSerialNumber', ...
                'PatientID', 'AccessionNumber', 'StudyID'};
            fields = [nameFields, wholeFields];
            tree = mlphibroker.DicomTree(root);
            names = {};
            for i = 1:height(tree.Series)
                h = tree.Series.FirstHeader{i};
                for f = fields
                    if ~isfield(h, f{1}) || isempty(h.(f{1})), continue; end
                    v = h.(f{1});
                    if isstruct(v)   % PN value representation (FamilyName, GivenName, ...)
                        parts = struct2cell(v); parts = parts(cellfun(@(x) ischar(x) || isstring(x), parts));
                        txt = strjoin(cellfun(@char, parts, 'UniformOutput', false), ' ');
                    else
                        txt = char(string(v));
                    end
                    txt = strtrim(txt);
                    if numel(txt) >= 3, names{end+1} = txt; end %#ok<AGROW>
                    if ismember(f{1}, nameFields)
                        words = regexp(txt, '[A-Za-z][A-Za-z''-]{2,}', 'match');
                        names = [names, words]; %#ok<AGROW>
                    end
                end
            end
            names = [names, cellstr(extra)];
            names = unique(strtrim(names(~cellfun(@isempty, names))), 'stable');
            % words that are PHI-free DICOM boilerplate and would only cause false quarantines
            stop = {'none', 'anonymous', 'unknown', 'brain', 'head', 'pet', 'dynamic', 'static', 'delay', ...
                'delayed', 'the', 'and', 'for', 'with', 'wholebody', 'whole', 'body'};
            names = names(~ismember(lower(names), stop));
            fid = fopen(outPath, 'w'); fprintf(fid, '%s\n', names{:}); fclose(fid);
            fprintf('deny list template: %d entries -> %s (PHI; review, add coordinator/investigator names)\n', numel(names), outPath);
        end

        function explainQuarantine(outRoot, denyListFile)
            %EXPLAINQUARANTINE  PHI-SIDE CONSOLE ONLY: for every quarantined sidecar under outRoot, print
            %   the key, the shape kind and the deny-list entry or matched text that caused the hit.
            %   Never paste this output anywhere outside the PHI zone; it shows the matched values.
            deny = {}; if nargin > 1 && ~isempty(denyListFile), deny = mlphibroker.SidecarPolicy.readDenyList(denyListFile); end
            p = mlphibroker.SidecarPolicy(deny);
            files = dir(fullfile(outRoot, '**', 'quarantine', '*.json'));
            if isempty(files), fprintf('no quarantined sidecars under %s\n', outRoot); return; end
            shapes = struct('date', '(?<![0-9])(19|20)\d{2}-?[01]\d-?[0-3]\d(?![0-9])|(?<![0-9])[01]?\d/[0-3]?\d/(19|20)\d{2}(?![0-9])', ...
                'clock_time', '(?<![0-9])[0-2]\d:[0-5]\d(:[0-5]\d(\.\d+)?)?(?![0-9])', 'dated_session', 'ses-\d{8,}');
            for i = 1:numel(files)
                fprintf('%s\n', files(i).name);
                sc = jsondecode(fileread(fullfile(files(i).folder, files(i).name)));
                for f = fieldnames(sc)'
                    txt = jsonencode(sc.(f{1}));
                    for k = fieldnames(shapes)'
                        m = regexp(txt, shapes.(k{1}), 'match');
                        if ~isempty(m), fprintf('   %-28s %-14s %s\n', f{1}, k{1}, strjoin(unique(m), ' ')); end
                    end
                    for d = deny
                        if ~isempty(regexpi(txt, ['(?<![A-Za-z])' regexptranslate('escape', d{1}) '(?![A-Za-z])'], 'once'))
                            fprintf('   %-28s %-14s entry "%s"\n', f{1}, 'deny_list', d{1});
                        end
                    end
                end
                h = p.scan(files(i).name); if h.any, fprintf('   %-28s file name carries a shape\n', ''); end
            end
        end

        function writeJson(path, s)
            fid = fopen(path, 'w'); fwrite(fid, jsonencode(s, 'PrettyPrint', true), 'char'); fclose(fid);
        end
    end
end
