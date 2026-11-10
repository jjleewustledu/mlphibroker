classdef PetSession
    %PETSESSION  One PET series -> Safe Harbor, idif2026-ready .nii.gz + .json.
    %   [n, dropped, quarantined] = mlphibroker.PetSession.convert(row, sub, ses, outDir, conv, policy, opts)
    %   Sidecar keys follow idif2026 loaders (TissueData: times/taus/timesMid; InputFuncData:
    %   halflife; decay_correction.stamp_decay_correction: DecayCorrected + DecayCorrection*) and
    %   PETUtilities.resolve_isotope (RadionuclideHalfLife, TracerRadionuclide, Radiopharmaceutical).
    %   The administration offset comes from the DICOM headers (dcm2niix omits
    %   RadiopharmaceuticalStartTime); every time is seconds from scan start, never a clock time.

    methods (Static)
        function [n, dropped, quarantined, paths] = convert(row, sub, ses, outDir, conv, policy, opts, run)
            import mlphibroker.*
            if nargin < 8 || isempty(run), run = 0; end
            n = 0; quarantined = 0; dropped = {}; paths = {};
            files = row.Files{1};
            admin = PetSession.adminOffset(files);
            prefix = sprintf('%s_%s_ser-', sub, ses);
            tmpOut = tempname(); mkdir(tmpOut);
            cleanup = onCleanup(@() rmdir(tmpOut, 's'));
            [jsons, niis] = conv.convertSeries(files, prefix, tmpOut);
            for k = 1:numel(jsons)
                src = jsondecode(fileread(jsons{k}));
                [kept, d] = policy.allow(src); dropped = union(dropped, d);
                [token, isotope, halflife_s] = TracerTable.resolve(kept, opts);
                [times, taus, timesMid] = PetSession.frameTiming(kept);
                sc = PetSession.sidecar(kept, sub, ses, token, isotope, halflife_s, times, taus, timesMid, admin);
                base = sprintf('%s_%s_trc-%s%s_proc-shdeid_pet', sub, ses, token, PetSession.runTag(run, k, numel(jsons)));
                hits = policy.scan([jsonencode(sc) ' ' base]);
                if hits.any
                    quarantined = quarantined + 1;
                    PetSession.quarantine(outDir, niis{k}, jsons{k}, base);
                    continue
                end
                Nifti1.scrubCopy(niis{k}, fullfile(outDir, [base '.nii.gz']));
                SidecarPolicy.writeJson(fullfile(outDir, [base '.json']), sc);
                paths{end+1} = fullfile(outDir, [base '.nii.gz']); %#ok<AGROW>
                n = n + 1;
            end
        end

        function admin = adminOffset(files)
            % Three states; header reads only; the clock values never leave this function.
            admin = struct('state', 'unmeasured', 'seconds', [], 'reason', '', 'weight_kg', []);
            acq = []; rst = [];
            for i = 1:numel(files)
                try
                    h = dicominfo(files{i}, 'UseDictionaryVR', true);
                catch
                    continue
                end
                if isempty(admin.weight_kg) && isfield(h, 'PatientWeight') && ~isempty(h.PatientWeight)
                    admin.weight_kg = double(h.PatientWeight);   % kg (DICOM 0010,1030)
                end
                if isfield(h, 'AcquisitionTime')
                    a = mlphibroker.PetSession.clockSeconds(h.AcquisitionTime);
                    if ~isempty(a), acq(end+1) = a; end %#ok<AGROW>
                end
                if isempty(rst) && isfield(h, 'RadiopharmaceuticalInformationSequence')
                    seq = h.RadiopharmaceuticalInformationSequence; fn = fieldnames(seq);
                    if ~isempty(fn) && isfield(seq.(fn{1}), 'RadiopharmaceuticalStartTime')
                        rst = mlphibroker.PetSession.clockSeconds(seq.(fn{1}).RadiopharmaceuticalStartTime);
                    end
                end
            end
            if isempty(acq), admin.reason = 'AcquisitionTime absent'; return; end  % weight may still be set
            if isempty(rst), admin.reason = 'RadiopharmaceuticalStartTime absent'; return; end
            d = rst - min(acq);
            if d > 43200, d = d - 86400; elseif d < -43200, d = d + 86400; end
            admin.state = 'measured'; admin.seconds = d;
        end

        function s = clockSeconds(t)
            t = strrep(strtrim(char(t)), ':', ''); s = [];
            tok = regexp(t, '^(\d{2})(\d{2})(\d{2}(\.\d+)?)?$', 'tokens', 'once');
            if isempty(tok), return; end
            sec = 0; if ~isempty(tok{3}), sec = str2double(tok{3}); end
            s = str2double(tok{1}) * 3600 + str2double(tok{2}) * 60 + sec;
        end

        function [times, taus, timesMid] = frameTiming(sc)
            if ~isfield(sc, 'FrameTimesStart') || ~isfield(sc, 'FrameDuration')
                error('mlphibroker:timing', 'FrameTimesStart and FrameDuration are required');
            end
            times = double(sc.FrameTimesStart(:))'; taus = double(sc.FrameDuration(:))';
            if numel(times) ~= numel(taus) || isempty(times)
                error('mlphibroker:timing', 'FrameTimesStart and FrameDuration lengths differ or are empty');
            end
            if any(taus <= 0) || times(1) < 0 || (numel(times) > 1 && any(diff(times) <= 0))
                error('mlphibroker:timing', 'frames must be positive in duration and strictly increasing');
            end
            if numel(times) > 1 && any(times(2:end) < times(1:end-1) + taus(1:end-1) - 1e-6)
                error('mlphibroker:timing', 'frames overlap');
            end
            timesMid = times + taus / 2;
        end

        function sc = sidecar(kept, sub, ses, token, isotope, halflife_s, times, taus, timesMid, admin)
            import mlphibroker.*
            sc = struct();
            sc.Subject = sub; sc.Session = ses; sc.Tracer = ['trc-' token]; sc.Modality = 'PET';
            sc.TracerRadionuclide = isotope; sc.RadionuclideHalfLife = halflife_s; sc.halflife = halflife_s;
            if isfield(kept, 'Radiopharmaceutical'), sc.Radiopharmaceutical = kept.Radiopharmaceutical;
            elseif isfield(kept, 'TracerName'), sc.Radiopharmaceutical = kept.TracerName; end   % idif2026 reads Radiopharmaceutical
            sc.times = times; sc.taus = taus; sc.timesMid = timesMid;
            sc.timeUnit = 'second'; sc.timeAnchor = 'scan_start';
            dc = ''; if isfield(kept, 'DecayCorrection'), dc = upper(char(kept.DecayCorrection)); end
            switch dc
                case 'START', zero_s = 0.0; corrected = true;
                case 'ADMIN'
                    if ~strcmp(admin.state, 'measured')
                        error('mlphibroker:decay', 'DecayCorrection is ADMIN but the administration offset is %s', admin.state);
                    end
                    zero_s = admin.seconds; corrected = true;
                case 'NONE', zero_s = 0.0; corrected = false;
                otherwise
                    error('mlphibroker:decay', 'DecayCorrection must be START, ADMIN or NONE; got "%s"', dc);
            end
            sc.DecayCorrected = corrected; sc.DecayCorrection = dc;
            sc.DecayCorrectionZeroTime = zero_s; sc.DecayCorrectionZeroTimeReference = 'scan_start';
            sc.DecayCorrectionZeroTimeUnits = 'second'; sc.DecayCorrectionHalflife = halflife_s;
            sc.DecayCorrectionRateConstant = log(2) / halflife_s; sc.DecayCorrectionRateConstantUnits = '1/second';
            sc.DecayCorrectionProvenance = [SafeHarbor.TOOL ' ' SafeHarbor.VERSION ' from DICOM DecayCorrection'];
            if ~isfield(kept, 'PatientWeight') && isfield(admin, 'weight_kg') && ~isempty(admin.weight_kg)
                sc.PatientWeight = admin.weight_kg;   % kg; released exact (operator ruling OQ-15 Option B, 2026-11-09)
                sc.PatientWeightUnits = 'kg';
            end
            sc.AdminFromScanStartState = admin.state;
            if strcmp(admin.state, 'measured'), sc.AdminFromScanStart_s = admin.seconds; else, sc.AdminFromScanStart_s = []; end
            pass = setdiff(SidecarPolicy.ALLOW, {'Radiopharmaceutical', 'RadionuclideHalfLife', 'FrameTimesStart', 'FrameDuration', 'Modality'});
            for i = 1:numel(pass)
                if isfield(kept, pass{i}), sc.(pass{i}) = kept.(pass{i}); end
            end
            sc.SafeHarbor = struct('tool', SafeHarbor.TOOL, 'version', SafeHarbor.VERSION, 'work_order', 'WO-0002r0', ...
                'policy', 'allow_list_then_relative_times', 'calendar_dates', 'none', 'clock_times', 'none', ...
                'nifti_header_scrubbed', true);
        end

        function tag = runTag(run, k, nOut)
            % run-NN names the study within a multi-study session; a series that dcm2niix splits
            % into several outputs appends the output index so names never collide.
            tag = '';
            if run > 0 && nOut > 1, tag = sprintf('_run-%02d%02d', run, k);
            elseif run > 0, tag = sprintf('_run-%02d', run);
            elseif nOut > 1, tag = sprintf('_run-%02d', k); end
        end

        function quarantine(outDir, nii, json, base)
            q = fullfile(outDir, 'quarantine'); if ~isfolder(q), mkdir(q); end
            copyfile(nii, fullfile(q, [base '.nii.gz'])); copyfile(json, fullfile(q, [base '.json']));
        end
    end
end
