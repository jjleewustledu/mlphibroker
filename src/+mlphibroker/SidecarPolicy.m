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
            'Radiopharmaceutical', 'TracerName', 'RadionuclideTotalDose', 'RadionuclideHalfLife', ...
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

        function writeJson(path, s)
            fid = fopen(path, 'w'); fwrite(fid, jsonencode(s, 'PrettyPrint', true), 'char'); fclose(fid);
        end
    end
end
