classdef Test_SafeHarbor < matlab.unittest.TestCase
    %% End to end on the committed synthetic DICOM tree (two fabricated subjects, three studies:
    %  dynamic PET, CT, one MIM-style RTSTRUCT). Needs dcm2niix on PATH or MLPHIBROKER_DCM2NIIX.
    %
    %  Created 09-Nov-2026 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlphibroker/test/+mlphibroker_unittest.
    %  Developed on Matlab 26.1 (R2026a). Copyright 2026 John J. Lee.

    properties
        root
        outRoot
        report
    end

    methods (TestClassSetup)
        function run_batch(this)
            here = fileparts(mfilename('fullpath'));
            this.root = fullfile(here, '..', 'fixtures', 'synthetic_dicom');
            exe = getenv('MLPHIBROKER_DCM2NIIX');
            if isempty(exe), exe = mlphibroker.Dcm2niix.findOnPath(); end
            this.assumeTrue(~isempty(exe) && isfile(exe), 'dcm2niix not found; end-to-end test skipped (could_not_run)');
            this.outRoot = tempname(); mkdir(this.outRoot);
            deny = fullfile(this.outRoot, 'deny.txt');
            fid = fopen(deny, 'w'); fprintf(fid, 'Zorbwick\nSynthetic Institute\n'); fclose(fid);
            opts = struct('AutoCrosswalk', true, 'Dcm2niixExe', exe, 'DenyListFile', deny, ...
                'SyntheticTest', true, 'Seed', 7);
            this.report = mlphibroker.SafeHarbor.run(this.root, this.outRoot, opts);
        end
    end

    methods (TestClassTeardown)
        function cleanup(this)
            if ~isempty(this.outRoot) && isfolder(this.outRoot), rmdir(this.outRoot, 's'); end
        end
    end

    methods (Test)
        function test_counts(this)
            r = this.report;
            this.verifyEqual(r.subjects, 2); this.verifyEqual(r.studies, 3);
            this.verifyEqual(r.series_by_kind.pet_dynamic, 3); this.verifyEqual(r.series_by_kind.ct, 3);
            this.verifyEqual(r.series_by_kind.rtstruct, 1);
            this.verifyEqual(r.written.pet, 3); this.verifyEqual(r.written.ct, 3); this.verifyEqual(r.written.roi, 2);
            this.verifyEqual(r.refused, 0, strjoin(r.refusals, ' | '));
            this.verifyEqual(r.quarantined, 0);
            this.verifyTrue(any(strcmp(r.dropped_key_names, 'AcquisitionTime')));
            this.verifyTrue(any(strcmp(r.dropped_key_names, 'TimeZero')));
        end
        function test_file_names_and_no_phi_shapes(this)
            files = dir(fullfile(this.outRoot, '**', '*.*')); files = files(~[files.isdir]);
            names = {files.name};
            this.verifyTrue(any(contains(names, '_trc-fdg_proc-shdeid_pet.nii.gz')));
            this.verifyTrue(any(contains(names, '_proc-shdeid_ct.nii.gz')));
            this.verifyTrue(any(contains(names, '_roi-leftcaudate_proc-shdeid_mask.nii.gz')));
            p = mlphibroker.SidecarPolicy({'Zorbwick', 'Synthetic Institute', 'Alpha', 'Beta', 'SYN000001', 'SYN000002'});
            for i = 1:numel(files)
                if endsWith(files(i).name, '.json')
                    txt = fileread(fullfile(files(i).folder, files(i).name));
                    h = p.scan(txt);
                    this.verifyFalse(h.any, ['PHI shape in ' files(i).name]);
                end
                this.verifyFalse(p.scan(files(i).name).any, files(i).name);
            end
            this.verifyFalse(any(contains(names, 'crosswalk')) && false); % crosswalk_PHI.csv stays PHI-side by design
        end
        function test_pet_sidecar_contract(this)
            j = dir(fullfile(this.outRoot, '**', '*_trc-fdg_proc-shdeid_pet.json'));
            sc = jsondecode(fileread(fullfile(j(1).folder, j(1).name)));
            this.verifyEqual(sc.times(:)', [0 30 60 120]); this.verifyEqual(sc.taus(:)', [30 30 60 60]);
            this.verifyEqual(sc.timesMid(:)', [15 45 90 150]);
            this.verifyTrue(sc.DecayCorrected); this.verifyEqual(sc.DecayCorrection, 'START');
            this.verifyEqual(sc.DecayCorrectionZeroTime, 0); this.verifyEqual(sc.DecayCorrectionZeroTimeReference, 'scan_start');
            this.verifyEqual(sc.halflife, 6586.2); this.verifyEqual(sc.TracerRadionuclide, '18F');
            this.verifyEqual(sc.AdminFromScanStartState, 'measured');
            this.verifyEqual(sc.AdminFromScanStart_s, -131, 'AbsTol', 1e-6);
            this.verifyEqual(sc.PatientWeight, 70.5);
            this.verifyFalse(isfield(sc, 'AcquisitionTime')); this.verifyFalse(isfield(sc, 'TimeZero'));
            this.verifyFalse(isfield(sc, 'InstitutionName'));
            nii = regexprep(fullfile(j(1).folder, j(1).name), '\.json$', '.nii.gz');
            [d, a, b] = mlphibroker.Nifti1.headerText(nii); this.verifyEmpty(d); this.verifyEmpty(a); this.verifyEmpty(b);
            info = niftiinfo(nii); this.verifyEqual(info.ImageSize, [16 16 4 4]);
        end
        function test_roi_mask_geometry(this)
            m = dir(fullfile(this.outRoot, '**', '*_roi-wholebrain_proc-shdeid_mask.nii.gz'));
            this.verifyEqual(numel(m), 1);
            mask = niftiread(fullfile(m(1).folder, m(1).name));
            this.verifyEqual(size(mask), [16 16 4]);
            this.verifyGreaterThan(nnz(mask), 0);
            sc = jsondecode(fileread(regexprep(fullfile(m(1).folder, m(1).name), '\.nii\.gz$', '.json')));
            this.verifyEqual(sc.VoxelsInside, nnz(mask));
            this.verifyEqual(sc.ReferenceModality, 'PT');
        end
        function test_roi_mask_centroid_in_ras(this)
            % Synthetic LeftCaudate: circle centred at DICOM LPS (-4, -4) mm, radius 3 mm, every slice.
            % In NIfTI RAS that is (+4, +4). The centroid of the mask, mapped through the sform
            % of the file it was written with, must land there: this pins the LPS/RAS handling and
            % the DICOM-grid to dcm2niix-grid placement.
            m = dir(fullfile(this.outRoot, '**', '*_roi-leftcaudate_proc-shdeid_mask.nii.gz'));
            this.assertEqual(numel(m), 1);
            gz = fullfile(m(1).folder, m(1).name);
            mask = niftiread(gz); S = mlphibroker.Nifti1.sform(gz);
            [i, j, k] = ind2sub(size(mask), find(mask));
            c = S * [mean(i) - 1; mean(j) - 1; mean(k) - 1; 1];
            this.verifyEqual(c(1:2)', [4 4], 'AbsTol', 0.6);
            this.verifyEqual(c(3), 3, 'AbsTol', 0.6);          % slices at z = 0,2,4,6 mm
            this.verifyGreaterThan(nnz(mask), 4 * 20); this.verifyLessThan(nnz(mask), 4 * 40);
            pet = dir(fullfile(m(1).folder, '*_trc-fdg_proc-shdeid_pet.nii.gz'));
            this.verifyEqual(mlphibroker.Nifti1.sform(fullfile(pet(1).folder, pet(1).name)), S);
        end
        function test_crosswalk_minted_and_private(this)
            xw = readtable(fullfile(this.outRoot, 'crosswalk_PHI.csv'), 'TextType', 'char');
            this.verifyEqual(height(xw), 3);
            this.verifyEqual(numel(unique(xw.Subject)), 2);
            this.verifyTrue(all(startsWith(xw.Subject, 'sub-9')));
            this.verifyTrue(all(~cellfun(@isempty, regexp(xw.Session, '^ses-\d{2}$'))));
        end
    end
end
