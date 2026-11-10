classdef PetSession
    %PETSESSION  One PET series -> Safe Harbor, idif2026-ready .nii.gz + .json.
    %   [n, dropped, quarantined] = mlphibroker.PetSession.convert(row, sub, ses, outDir, conv, policy, opts)
    %   Sidecar keys follow idif2026 loaders (TissueData: times/taus/timesMid; InputFuncData:
    %   halflife; decay_correction.stamp_decay_correction: DecayCorrected + DecayCorrection*) and
    %   PETUtilities.resolve_isotope (RadionuclideHalfLife, TracerRadionuclide, Radiopharmaceutical).
    %   The administration offset comes from the DICOM headers (dcm2niix omits
    %   RadiopharmaceuticalStartTime); every time is seconds from scan start, never a clock time.

    methods (Static)
        function [n, dropped, quarantined, paths, detail] = convert(row, sub, ses, outDir, conv, policy, opts, run)
            import mlphibroker.*
            if nargin < 8 || isempty(run), run = 0; end
            n = 0; quarantined = 0; dropped = {}; paths = {}; detail = {};
            files = row.Files{1};
            admin = PetSession.adminOffset(files);
            prefix = sprintf('%s_%s_ser-', sub, ses);
            tmpOut = tempname(); mkdir(tmpOut);
            cleanup = onCleanup(@() rmdir(tmpOut, 's'));
            [jsons, niis] = conv.convertSeries(files, prefix, tmpOut);
            for k = 1:numel(jsons)
                src = jsondecode(fileread(jsons{k}));
                [kept, d] = policy.allow(src); dropped = union(dropped, d);
                % dcm2niix (esp. <= 2023, and for single-frame PET) omits keys the DICOM headers hold:
                if (~isfield(kept, 'RadionuclideHalfLife') || isempty(kept.RadionuclideHalfLife)) && ~isempty(admin.halflife_s)
                    kept.RadionuclideHalfLife = admin.halflife_s; kept.RadionuclideHalfLifeSource = 'dicom_headers';
                end
                if ~isfield(kept, 'Radiopharmaceutical') && ~isfield(kept, 'TracerName') && ~isempty(admin.radiopharmaceutical)
                    kept.Radiopharmaceutical = admin.radiopharmaceutical;
                end
                [token, isotope, halflife_s] = TracerTable.resolve(kept, opts);
                nuclideCheck = 'absent';
                if isfield(kept, 'TracerRadionuclide') && ~isempty(kept.TracerRadionuclide)
                    % dcm2niix copies the DICOM radionuclide code meaning, e.g. "^68^Gallium", "Gallium-68",
                    % "Ga-68" or "68Ga": compare mass number and element, never the raw string.
                    [okN, nuclideCheck] = PetSession.nuclideAgrees(kept.TracerRadionuclide, isotope);
                    if ~okN
                        error('mlphibroker:isotope', 'dcm2niix TracerRadionuclide (%s) disagrees with the isotope selected by the half-life (%s)', nuclideCheck, isotope);
                    end
                end
                if isfield(kept, 'FrameTimesStart') && isfield(kept, 'FrameDuration')
                    [times, taus, timesMid] = PetSession.frameTiming(kept); timingSource = 'dcm2niix_sidecar';
                else
                    [times, taus, timesMid, timingSource] = PetSession.frameTimingFromHeaders(admin.frames);
                end
                sc = PetSession.sidecar(kept, sub, ses, token, isotope, halflife_s, times, taus, timesMid, admin);
                sc.FrameTimingSource = timingSource;
                sc.TracerRadionuclideCheck = nuclideCheck;   % 'agrees' | 'unparsed:<text>' | 'absent'
                base = sprintf('%s_%s_trc-%s%s_proc-shdeid_pet', sub, ses, token, PetSession.runTag(run, k, numel(jsons)));
                [hits, byField] = policy.scanStruct(sc, base);
                if hits.any
                    quarantined = quarantined + 1;
                    PetSession.quarantine(outDir, niis{k}, jsons{k}, base);
                    detail{end+1} = sprintf('quarantined: %s', strjoin(byField, ', ')); %#ok<AGROW>
                    continue
                end
                Nifti1.scrubCopy(niis{k}, fullfile(outDir, [base '.nii.gz']));
                SidecarPolicy.writeJson(fullfile(outDir, [base '.json']), sc);
                paths{end+1} = fullfile(outDir, [base '.nii.gz']); %#ok<AGROW>
                detail{end+1} = sprintf('written: %d frames, timing %s, admin %s', numel(times), timingSource, admin.state); %#ok<AGROW>
                n = n + 1;
            end
        end

        function admin = adminOffset(files)
            %ADMINOFFSET  Header facts for one PET series (DICOM headers only; clock values never leave).
            %   admin.state/seconds/reason : administration time relative to scan start (three states)
            %   admin.weight_kg            : PatientWeight when present
            %   admin.halflife_s, admin.radiopharmaceutical : from RadiopharmaceuticalInformationSequence
            %   admin.frames               : per-frame table (acq_s, duration_s, ref_ms) for frame timing
            admin = struct('state', 'unmeasured', 'seconds', [], 'reason', '', 'weight_kg', [], ...
                'halflife_s', [], 'radiopharmaceutical', '', 'frames', []);
            acq = []; rst = []; acqNo = []; dur = []; ref = [];
            for i = 1:numel(files)
                try
                    h = dicominfo(files{i}, 'UseDictionaryVR', true);
                catch
                    continue
                end
                if isempty(admin.weight_kg) && isfield(h, 'PatientWeight') && ~isempty(h.PatientWeight)
                    admin.weight_kg = double(h.PatientWeight);   % kg (DICOM 0010,1030)
                end
                a = NaN;
                if isfield(h, 'AcquisitionTime')
                    a0 = mlphibroker.PetSession.clockSeconds(h.AcquisitionTime);
                    if ~isempty(a0), a = a0; acq(end+1) = a0; end %#ok<AGROW>
                end
                if isfield(h, 'RadiopharmaceuticalInformationSequence')
                    seq = h.RadiopharmaceuticalInformationSequence; fn = fieldnames(seq);
                    if ~isempty(fn)
                        it = seq.(fn{1});
                        if isempty(rst) && isfield(it, 'RadiopharmaceuticalStartTime')
                            rst = mlphibroker.PetSession.clockSeconds(it.RadiopharmaceuticalStartTime);
                        end
                        if isempty(admin.halflife_s) && isfield(it, 'RadionuclideHalfLife') && ~isempty(it.RadionuclideHalfLife)
                            admin.halflife_s = double(it.RadionuclideHalfLife);
                        end
                        if isempty(admin.radiopharmaceutical) && isfield(it, 'Radiopharmaceutical') && ~isempty(it.Radiopharmaceutical)
                            admin.radiopharmaceutical = char(string(it.Radiopharmaceutical));
                        end
                    end
                end
                n = NaN; d = NaN; r = NaN;
                if isfield(h, 'AcquisitionNumber') && ~isempty(h.AcquisitionNumber), n = double(h.AcquisitionNumber); end
                if isfield(h, 'ActualFrameDuration') && ~isempty(h.ActualFrameDuration), d = double(h.ActualFrameDuration) / 1000; end
                if isfield(h, 'FrameReferenceTime') && ~isempty(h.FrameReferenceTime), r = double(h.FrameReferenceTime); end
                acqNo(end+1) = n; dur(end+1) = d; ref(end+1) = r; %#ok<AGROW>
                if numel(acqNo) > numel(acq), acq(end+1) = NaN; end %#ok<AGROW>
            end
            if ~isempty(acqNo)
                admin.frames = table(acqNo', acq(1:numel(acqNo))', dur', ref', 'VariableNames', {'acq_no', 'acq_s', 'duration_s', 'ref_ms'});
            end
            good = acq(~isnan(acq));
            if isempty(good), admin.reason = 'AcquisitionTime absent'; return; end
            if isempty(rst), admin.reason = 'RadiopharmaceuticalStartTime absent'; return; end
            d = rst - min(good);
            if d > 43200, d = d - 86400; elseif d < -43200, d = d + 86400; end
            admin.state = 'measured'; admin.seconds = d;
        end

        function [times, taus, timesMid, source] = frameTimingFromHeaders(fr)
            %FRAMETIMINGFROMHEADERS  Frame starts/durations (s) from per-file headers, grouped by frame.
            %   Frames are grouped by AcquisitionNumber (else FrameReferenceTime); start = acquisition
            %   time of the frame minus that of the first frame; duration = ActualFrameDuration.
            %   A single-frame (static or delayed) scan gives times = 0, taus = its duration.
            source = 'dicom_headers';
            if isempty(fr), error('mlphibroker:timing', 'no DICOM frame headers'); end
            key = fr.acq_no; if all(isnan(key)), key = fr.ref_ms; end
            if all(isnan(key)), key = ones(height(fr), 1); end
            [keys, ~, g] = unique(key(~isnan(key)));
            frv = fr(~isnan(key), :);
            times = zeros(1, numel(keys)); taus = zeros(1, numel(keys)); refs = nan(1, numel(keys));
            for k = 1:numel(keys)
                rows = frv(g == k, :);
                a = rows.acq_s(~isnan(rows.acq_s)); d = rows.duration_s(~isnan(rows.duration_s)); r = rows.ref_ms(~isnan(rows.ref_ms));
                if isempty(a) || isempty(d), error('mlphibroker:timing', 'frame %d lacks AcquisitionTime or ActualFrameDuration', k); end
                times(k) = min(a); taus(k) = mode(d); if ~isempty(r), refs(k) = mode(r); end
            end
            [times, order] = sort(times); taus = taus(order); refs = refs(order);
            times = times - times(1);
            times(times < -43200) = times(times < -43200) + 86400;     % midnight wrap
            if any(~isnan(refs))   % FrameReferenceTime is the frame midpoint in ms: cross-check
                mid = times + taus / 2;
                ok = isnan(refs) | abs(refs / 1000 - mid) <= max(1, 0.05 * taus);
                if ~all(ok)
                    error('mlphibroker:timing', 'FrameReferenceTime disagrees with AcquisitionTime/ActualFrameDuration on %d frame(s)', sum(~ok));
                end
            end
            if any(taus <= 0) || (numel(times) > 1 && any(diff(times) <= 0))
                error('mlphibroker:timing', 'frames must be positive in duration and strictly increasing');
            end
            timesMid = times + taus / 2;
        end

        function [ok, how] = nuclideAgrees(txt, isotope)
            %NUCLIDEAGREES  Does a radionuclide string ("^68^Gallium", "Ga-68", "68Ga") name `isotope` ("68Ga")?
            %   ok is true when the mass numbers agree and the element (symbol or English name) agrees,
            %   and also when the string cannot be parsed (then how = 'unparsed:<text>'): the half-life
            %   already fixed the isotope, and an unreadable code meaning is recorded, not fatal.
            elements = { 'ga', 'gallium'; 'cu', 'copper'; 'f', 'fluorine'; 'o', 'oxygen'; 'c', 'carbon'; ...
                'zr', 'zirconium'; 'rb', 'rubidium'; 'n', 'nitrogen'; 'i', 'iodine'; 'ge', 'germanium'; ...
                'na', 'sodium'; 'cs', 'caesium'; 'cs', 'cesium'; 'tc', 'technetium'; 'y', 'yttrium'; 'lu', 'lutetium' };
            t = lower(char(string(txt)));
            massT = regexp(t, '\d+', 'match', 'once');
            elT = regexp(t, '[a-z]+', 'match', 'once');
            iso = lower(char(isotope));
            massI = regexp(iso, '\d+', 'match', 'once'); symI = regexp(iso, '[a-z]+', 'match', 'once');
            if isempty(massT) || isempty(elT)
                ok = true; how = ['unparsed:' regexprep(char(string(txt)), '[^A-Za-z0-9^-]', '')]; return
            end
            r = find(strcmp(elements(:, 1), symI));
            names = [{symI}, elements(r, 2)'];
            ok = strcmp(massT, massI) && any(strcmp(elT, names));
            if ok, how = 'agrees'; else, how = regexprep(char(string(txt)), '[^A-Za-z0-9^-]', ''); end
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
            if isfield(kept, 'RadionuclideHalfLifeSource'), sc.RadionuclideHalfLifeSource = kept.RadionuclideHalfLifeSource; end
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
