classdef Test_Nifti1 < matlab.unittest.TestCase
    %% Header scrub zeroes descrip/aux_file/db_name and leaves voxel bytes untouched.
    %
    %  Created 09-Nov-2026 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlphibroker/test/+mlphibroker_unittest.
    %  Developed on Matlab 26.1 (R2026a). Copyright 2026 John J. Lee.

    methods (Test)
        function test_scrub_preserves_data(this)
            tmp = tempname(); mkdir(tmp);
            V = int16(reshape(1:(4*5*6*3), 4, 5, 6, 3));
            nii = fullfile(tmp, 'v.nii');
            niftiwrite(V, nii);
            % plant a clock time in descrip as dcm2niix does
            fid = fopen(nii, 'r+'); fseek(fid, 148, 'bof'); fwrite(fid, 'Time=120011.000', 'char'); fclose(fid);
            gz = gzip(nii, tmp); gz = gz{1};
            [d0, ~, ~] = mlphibroker.Nifti1.headerText(gz); this.verifyEqual(d0, 'Time=120011.000');
            out = fullfile(tmp, 'scrubbed.nii.gz');
            mlphibroker.Nifti1.scrubCopy(gz, out);
            [d1, a1, b1] = mlphibroker.Nifti1.headerText(out);
            this.verifyEmpty(d1); this.verifyEmpty(a1); this.verifyEmpty(b1);
            this.verifyEqual(niftiread(out), V);
            info = niftiinfo(out); this.verifyEqual(info.ImageSize, size(V));
            rmdir(tmp, 's');
        end
        function test_refuses_non_nifti(this)
            tmp = tempname(); mkdir(tmp);
            f = fullfile(tmp, 'x.bin'); fid = fopen(f, 'w'); fwrite(fid, zeros(1, 400, 'uint8')); fclose(fid);
            gz = gzip(f, tmp); gz = gz{1};
            this.verifyError(@() mlphibroker.Nifti1.scrubCopy(gz, fullfile(tmp, 'y.nii.gz')), 'mlphibroker:nifti');
            rmdir(tmp, 's');
        end
    end
end
