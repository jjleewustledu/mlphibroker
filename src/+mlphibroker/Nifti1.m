classdef Nifti1
    %NIFTI1  Byte-level NIfTI-1 .nii.gz helpers: header scrub, sform read, write-like-reference.
    %   mlphibroker.Nifti1.scrubCopy(inGz, outGz) zeroes db_name (offset 14, 18 B), descrip
    %   (148, 80 B) and aux_file (228, 24 B); dcm2niix writes the acquisition clock time into
    %   descrip as "Time=hhmmss". Refuses anything but a little-endian NIfTI-1 single file.
    %   S = mlphibroker.Nifti1.sform(gz) returns the 4x4 voxel(0-based)->RAS matrix (srow_*).
    %   mlphibroker.Nifti1.writeLike(refGz, outGz, V) writes uint8 V (3-D) with the reference
    %   file's header (same grid and sform), datatype uint8, scrubbed text fields.
    %   Uses MATLAB's native gunzip/gzip (Java byte arrays are passed by value in MATLAB, so a
    %   streaming Java gzip read silently yields zeros; do not reintroduce it).

    methods (Static)
        function scrubCopy(inGz, outGz)
            tmp = tempname(); mkdir(tmp);
            cleanup = onCleanup(@() rmdir(tmp, 's'));
            raw = gunzip(inGz, tmp); raw = raw{1};
            mlphibroker.Nifti1.assertNifti1(raw);
            fid = fopen(raw, 'r+');
            for range = {[14, 18], [148, 80], [228, 24]}
                fseek(fid, range{1}(1), 'bof'); fwrite(fid, zeros(1, range{1}(2), 'uint8'), 'uint8');
            end
            fclose(fid);
            gz = gzip(raw, tmp); gz = gz{1};
            movefile(gz, outGz, 'f');
        end

        function [descrip, auxfile, dbname] = headerText(gzPath)
            b = mlphibroker.Nifti1.headerBytes(gzPath);
            descrip = deblank(char(b(149:228))'); auxfile = deblank(char(b(229:252))'); dbname = deblank(char(b(15:32))');
        end

        function S = sform(gzPath)
            b = mlphibroker.Nifti1.headerBytes(gzPath);
            code = typecast(b(255:256), 'int16');
            if code <= 0, error('mlphibroker:nifti', 'reference NIfTI has no sform'); end
            rows = typecast(b(281:328), 'single');     % srow_x, srow_y, srow_z (12 floats)
            S = [double(reshape(rows, 4, 3))'; 0 0 0 1];
        end

        function dims = dim(gzPath)
            b = mlphibroker.Nifti1.headerBytes(gzPath);
            d = typecast(b(41:56), 'int16');
            dims = double(d(2:1 + d(1)))';
        end

        function writeLike(refGz, outGz, V)
            % V: 3-D uint8 array on the reference grid (dims 1-3 of the reference).
            b = mlphibroker.Nifti1.headerBytes(refGz);
            d = typecast(b(41:56), 'int16');
            if ~isequal(double(d(2:4))', size(V, 1:3))
                error('mlphibroker:nifti', 'mask size %s does not match the reference grid %s', ...
                    mat2str(size(V)), mat2str(double(d(2:4))'));
            end
            d(1) = 3; d(5:8) = 1;
            b(41:56) = typecast(int16(d), 'uint8');
            b(71:72) = typecast(int16(2), 'uint8');       % datatype uint8
            b(73:74) = typecast(int16(8), 'uint8');       % bitpix
            b(69:70) = typecast(int16(0), 'uint8');       % intent_code
            b(109:112) = typecast(single(352), 'uint8');  % vox_offset
            b(113:116) = typecast(single(1), 'uint8');    % scl_slope
            b(117:120) = typecast(single(0), 'uint8');    % scl_inter
            b(125:132) = 0;                               % cal_max, cal_min
            b(15:32) = 0; b(149:228) = 0; b(229:252) = 0;  % db_name, descrip, aux_file
            b(349:352) = 0;                               % extension flag
            tmp = tempname(); mkdir(tmp);
            cleanup = onCleanup(@() rmdir(tmp, 's'));
            raw = fullfile(tmp, 'mask.nii');
            fid = fopen(raw, 'w'); fwrite(fid, b, 'uint8'); fwrite(fid, uint8(V(:)), 'uint8'); fclose(fid);
            gz = gzip(raw, tmp); gz = gz{1};
            movefile(gz, outGz, 'f');
        end
    end

    methods (Static, Access = private)
        function b = headerBytes(gzPath)
            tmp = tempname(); mkdir(tmp);
            cleanup = onCleanup(@() rmdir(tmp, 's'));
            raw = gunzip(gzPath, tmp); raw = raw{1};
            mlphibroker.Nifti1.assertNifti1(raw);
            fid = fopen(raw, 'r'); b = fread(fid, 352, 'uint8=>uint8'); fclose(fid);
        end

        function assertNifti1(rawPath)
            fid = fopen(rawPath, 'r');
            if fid < 0, error('mlphibroker:nifti', 'cannot open the uncompressed NIfTI'); end
            hdr = fread(fid, 1, 'int32=>double');
            fseek(fid, 344, 'bof'); magic = fread(fid, 3, 'uint8=>char')';
            fclose(fid);
            if hdr ~= 348 || ~strcmp(magic, 'n+1')
                error('mlphibroker:nifti', 'not a little-endian NIfTI-1 single file; refusing');
            end
        end
    end
end
