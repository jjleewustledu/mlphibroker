classdef RtStruct
    %RTSTRUCT  MIM RTSTRUCT -> one Safe Harbor NIfTI mask per ROI, on the converted image grid.
    %   [n, quarantined] = mlphibroker.RtStruct.convert(rsRow, studyRows, written, sub, ses, outDir, policy)
    %   written: containers.Map SeriesInstanceUID -> path of that series' converted .nii.gz.
    %   The reference series is the one named in RTReferencedSeriesSequence, else the first PET or
    %   CT sharing the FrameOfReferenceUID. Contours are rasterised with dicomContours/createMask
    %   (Image Processing Toolbox) on the reference's spatial referencing from dicomreadVolume,
    %   using the FIRST TIME FRAME only for a dynamic PET; the mask is then placed voxel by voxel
    %   onto the dcm2niix output grid through that file's sform (DICOM LPS <-> NIfTI RAS), and
    %   written with the reference file's own header (Nifti1.writeLike), so the mask and the PET
    %   the kinetics use share one grid and one affine. ROI names pass the deny-list and PHI-shape
    %   scan, then become roi-<name> tokens (alphanumerics only). No dates, times, names or UIDs.

    methods (Static)
        function [n, quarantined] = convert(rsRow, studyRows, written, sub, ses, outDir, policy, runOf)
            import mlphibroker.*
            n = 0; quarantined = 0;
            info = dicominfo(rsRow.Files{1}{1}, 'UseDictionaryVR', true);
            ref = RtStruct.referenceSeries(info, rsRow, studyRows);
            if isempty(ref), error('mlphibroker:rtstruct', 'no referenced PET or CT series in this study'); end
            uid = ref.SeriesInstanceUID{1};
            if ~isKey(written, uid), error('mlphibroker:rtstruct', 'referenced series was not converted'); end
            refGz = written(uid);
            run = 0; if nargin >= 8 && isa(runOf, 'containers.Map') && isKey(runOf, uid), run = runOf(uid); end
            files = RtStruct.firstFrameFiles(ref.Files{1});
            spatial = RtStruct.geometryFromHeaders(files);   % own geometry: slice order along the normal, any direction
            S = Nifti1.sform(refGz);
            dims = Nifti1.dim(refGz);
            contours = dicomContours(info);
            for i = 1:height(contours.ROIs)
                name = char(contours.ROIs.Name(i));
                token = lower(regexprep(name, '[^A-Za-z0-9]', ''));
                if isempty(token), token = sprintf('roi%d', i); end
                base = sprintf('%s_%s_roi-%s%s_proc-shdeid_mask', sub, ses, token, PetSession.runTag(run, 1, 1));
                sc = struct('Subject', sub, 'Session', ses, 'Modality', 'ROI', 'ROIName', name, ...
                    'ROINumber', double(contours.ROIs.Number(i)), 'ReferenceModality', char(ref.Modality), ...
                    'ReferenceSeriesNumber', ref.SeriesNumber, 'ReferenceFile', regexprep(refGz, '^.*[\\/]', ''), ...
                    'Source', 'RTSTRUCT (MIM)', 'MaskEncoding', 'uint8 1 inside', ...
                    'SafeHarbor', struct('tool', SafeHarbor.TOOL, 'version', SafeHarbor.VERSION, ...
                    'work_order', 'WO-0002r0', 'policy', 'allow_list', 'calendar_dates', 'none', ...
                    'clock_times', 'none', 'nifti_header_scrubbed', true));
                hits = policy.scan([jsonencode(sc) ' ' base]);
                if hits.any
                    quarantined = quarantined + 1;
                    continue
                end
                mask = RtStruct.rasterize(contours, i, spatial);
                onGrid = RtStruct.toReferenceGrid(logical(mask), spatial, S, dims(1:3));
                Nifti1.writeLike(refGz, fullfile(outDir, [base '.nii.gz']), uint8(onGrid));
                sc.VoxelsInside = nnz(onGrid); sc.VoxelsInsideDicomGrid = nnz(mask);
                SidecarPolicy.writeJson(fullfile(outDir, [base '.json']), sc);
                n = n + 1;
            end
        end

        function ref = referenceSeries(info, rsRow, studyRows)
            ref = [];
            wanted = '';
            try
                rf = info.ReferencedFrameOfReferenceSequence; f1 = fieldnames(rf);
                st = rf.(f1{1}).RTReferencedStudySequence; f2 = fieldnames(st);
                se = st.(f2{1}).RTReferencedSeriesSequence; f3 = fieldnames(se);
                wanted = se.(f3{1}).SeriesInstanceUID;
            catch
            end
            img = studyRows(ismember(studyRows.Kind, {'pet_dynamic', 'pet_static', 'ct'}), :);
            if isempty(img), return; end
            hit = strcmp(img.SeriesInstanceUID, wanted);
            if ~any(hit) && ~isempty(rsRow.FrameOfReferenceUID{1})
                hit = strcmp(img.FrameOfReferenceUID, rsRow.FrameOfReferenceUID{1});
            end
            if any(hit), ref = img(find(hit, 1), :); else, ref = img(1, :); end
        end

        function spatial = geometryFromHeaders(files)
            %GEOMETRYFROMHEADERS  Slice geometry of one 3-D volume from its DICOM headers.
            %   Returns PatientPositions [nS x 3] (LPS, mm) sorted ascending along the slice normal,
            %   PatientOrientations [2 x 3 x nS], PixelSpacings [nS x 2] ([row col] mm) and ImageSize
            %   [rows cols nS]. Works whatever order the files were stored in.
            files = cellstr(files);
            P = []; O = []; PS = []; RC = [];
            for i = 1:numel(files)
                h = dicominfo(files{i}, 'UseDictionaryVR', true);
                if ~isfield(h, 'ImagePositionPatient') || ~isfield(h, 'ImageOrientationPatient'), continue; end
                P(end+1, :) = double(h.ImagePositionPatient(:))'; %#ok<AGROW>
                o = double(h.ImageOrientationPatient(:))'; O(:, :, end+1) = [o(1:3); o(4:6)]; %#ok<AGROW>
                PS(end+1, :) = double(h.PixelSpacing(:))'; %#ok<AGROW>
                RC(end+1, :) = double([h.Rows, h.Columns]); %#ok<AGROW>
            end
            if isempty(P), error('mlphibroker:rtstruct', 'reference series has no ImagePositionPatient'); end
            O = O(:, :, 2:end);   % drop the empty first page created by end+1 on an empty array
            X = O(1, :, 1) / norm(O(1, :, 1)); Y = O(2, :, 1) / norm(O(2, :, 1)); N = cross(X, Y);
            [~, order] = sort(P * N');
            spatial = struct('PatientPositions', P(order, :), 'PatientOrientations', O(:, :, order), ...
                'PixelSpacings', PS(order, :), 'ImageSize', [RC(1, 1), RC(1, 2), numel(order)]);
        end

        function mask = rasterize(contours, roiIndex, spatial)
            %RASTERIZE  CLOSED_PLANAR contours of one ROI -> logical [rows cols slices] on `spatial`.
            %   Each contour's points (LPS mm) are projected onto the slice whose plane they lie in and
            %   filled with poly2mask; several contours on one slice combine by XOR (holes). Points
            %   off every slice plane by more than half a slice spacing refuse the ROI.
            P1 = spatial.PatientPositions(1, :);
            O = spatial.PatientOrientations(:, :, 1);
            X = O(1, :) / norm(O(1, :)); Y = O(2, :) / norm(O(2, :)); N = cross(X, Y);
            ps = spatial.PixelSpacings(1, :);
            nS = spatial.ImageSize(3);
            if nS > 1, dz = dot(spatial.PatientPositions(2, :) - P1, N); else, dz = 1; end
            mask = false(spatial.ImageSize);
            data = contours.ROIs.ContourData{roiIndex};
            if ~iscell(data), data = {data}; end
            for c = 1:numel(data)
                pts = data{c};
                if isempty(pts), continue; end
                if size(pts, 2) ~= 3, pts = reshape(pts, 3, [])'; end
                rel = pts - P1;
                zc = rel * N' / dz;                               % continuous slice index (0-based)
                sIdx = round(mean(zc)) + 1;
                if abs(mean(zc) - (sIdx - 1)) > 0.5 || sIdx < 1 || sIdx > nS
                    error('mlphibroker:rtstruct', 'a contour lies off every slice plane of the reference series');
                end
                col = rel * X' / ps(2) + 1;                       % 1-based pixel centre coordinates
                row = rel * Y' / ps(1) + 1;
                plane = poly2mask(col, row, spatial.ImageSize(1), spatial.ImageSize(2));
                mask(:, :, sIdx) = xor(mask(:, :, sIdx), plane);
            end
        end

        function files = firstFrameFiles(files)
            % Dynamic PET: keep the files of the earliest acquisition (one spatial volume).
            files = cellstr(files);
            acq = nan(1, numel(files));
            for i = 1:numel(files)
                try
                    h = dicominfo(files{i}, 'UseDictionaryVR', true);
                    if isfield(h, 'AcquisitionNumber'), acq(i) = double(h.AcquisitionNumber); end
                catch
                end
            end
            if all(isnan(acq)) || numel(unique(acq(~isnan(acq)))) <= 1, return; end
            files = files(acq == min(acq));
        end

        function onGrid = toReferenceGrid(mask, spatial, S, dims)
            % mask: [rows cols slices] logical from createMask (DICOM order, 1-based).
            % For each NIfTI voxel (0-based i,j,k) compute RAS = S*[i j k 1]', LPS = [-x -y z],
            % then its (row, col, slice) in the DICOM grid and copy the mask value (nearest).
            P1 = spatial.PatientPositions(1, :);
            O = spatial.PatientOrientations(:, :, 1);          % 2x3: X (along columns), Y (along rows)
            X = O(1, :) / norm(O(1, :)); Y = O(2, :) / norm(O(2, :)); N = cross(X, Y);
            ps = spatial.PixelSpacings(1, :);                  % [row spacing, col spacing]
            nS = size(spatial.PatientPositions, 1);
            if nS > 1
                dz = dot(spatial.PatientPositions(2, :) - P1, N);
            else
                dz = 1;
            end
            [I, J, K] = ndgrid(0:dims(1) - 1, 0:dims(2) - 1, 0:dims(3) - 1);
            ras = S * [I(:)'; J(:)'; K(:)'; ones(1, numel(I))];
            lps = [-ras(1, :); -ras(2, :); ras(3, :)];
            rel = lps - P1';
            c = round(dot(rel, repmat(X', 1, size(rel, 2))) / ps(2)) + 1;
            r = round(dot(rel, repmat(Y', 1, size(rel, 2))) / ps(1)) + 1;
            s = round(dot(rel, repmat(N', 1, size(rel, 2))) / dz) + 1;
            ok = r >= 1 & r <= size(mask, 1) & c >= 1 & c <= size(mask, 2) & s >= 1 & s <= size(mask, 3);
            onGrid = false(numel(I), 1);
            idx = sub2ind(size(mask), r(ok), c(ok), s(ok));
            onGrid(ok) = mask(idx);
            onGrid = reshape(onGrid, dims(1), dims(2), dims(3));
        end
    end
end
