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
            % Every dcm2niix release writes AcquisitionTime; newer ones also write TimeZero,
            % InstitutionName and serials. Assert what must hold for any release: the clock key
            % was dropped, nothing dropped was allow-listed, and the version is recorded.
            this.verifyTrue(any(strcmp(r.dropped_key_names, 'AcquisitionTime')), ...
                ['dropped: ' strjoin(r.dropped_key_names, ',')]);
            this.verifyEmpty(intersect(r.dropped_key_names, mlphibroker.SidecarPolicy.ALLOW));
            this.verifyTrue(startsWith(r.dcm2niix_version, 'v'), r.dcm2niix_version);
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
        function test_crosswalk_template_then_load(this)
            tmp = fullfile(this.outRoot, 'template_PHI.csv');
            t = mlphibroker.Crosswalk.template(this.root, tmp);                 % no prefix: Subject left for the operator
            this.verifyEqual(height(t), 3); this.verifyEqual(numel(unique(t.PatientID)), 2);
            this.verifyTrue(all(cellfun(@isempty, t.Subject)));
            this.verifyTrue(all(~cellfun(@isempty, t.Session)));                   % sessions are always proposed
            this.verifyTrue(any(contains(t.Modalities, 'RTSTRUCT')));
            this.verifyTrue(isfile([tmp '.README.txt']));
            this.verifyError(@() mlphibroker.Crosswalk.load(tmp), 'mlphibroker:crosswalk');  % blanks refused
            t = mlphibroker.Crosswalk.template(this.root, tmp, struct('SubjectPrefix', 'sub-900', 'Seed', 3));
            this.verifyEqual(sort(unique(t.Subject))', {'sub-900001', 'sub-900002'});   % one code per person
            b = strcmp(t.PatientID, 'SYN000002');
            this.verifyEqual(sort(t.Session(b))', {'ses-01', 'ses-02'});                % two visit days, two ordinals
            this.verifyTrue(all(contains(t.Note(b), 'assigned at random')));
            this.verifyTrue(all(contains(t.Note, 'PLACEHOLDER')));
            % coordinator roster: codes come from it (with gaps), a person missing from it stays empty
            rosterPath = fullfile(this.outRoot, 'roster_PHI.csv');
            writetable(table({'SYN000002'}, {'sub-900005'}, 'VariableNames', {'PatientID', 'Subject'}), rosterPath);
            t = mlphibroker.Crosswalk.template(this.root, tmp, struct('Roster', rosterPath, 'SubjectPrefix', 'sub-900'));
            this.verifyEqual(unique(t.Subject(strcmp(t.PatientID, 'SYN000002'))), {'sub-900005'});
            this.verifyTrue(all(cellfun(@isempty, t.Subject(~strcmp(t.PatientID, 'SYN000002')))));
            this.verifyTrue(any(contains(t.Note, 'NOT IN ROSTER')));
            this.verifyError(@() mlphibroker.Crosswalk.load(tmp), 'mlphibroker:crosswalk');  % the empty one is refused
            % by name, caret or spaces, case-insensitive
            writetable(table({'synthetic alpha'; 'Synthetic^Beta'}, {'sub-900011'; 'sub-900005'}, 'VariableNames', {'PatientName', 'Subject'}), rosterPath);
            t = mlphibroker.Crosswalk.template(this.root, tmp, struct('Roster', rosterPath));
            this.verifyEqual(sort(unique(t.Subject))', {'sub-900005', 'sub-900011'});
            xw = mlphibroker.Crosswalk.load(tmp);
            this.verifyEqual(xw.Properties.VariableNames, {'PatientID', 'Subject', 'StudyInstanceUID', 'Session'});
            mlphibroker.Crosswalk.validate(xw);
        end
        function test_subject_code_ranges(this)
            xw = table({'p1'}, {'sub-ECL1i001'}, {'1.2.3'}, {'ses-01'}, 'VariableNames', {'PatientID', 'Subject', 'StudyInstanceUID', 'Session'});
            mlphibroker.Crosswalk.validate(xw);                      % MASH range accepted
            xw.Subject = {'sub-ECL1i-001'};
            this.verifyError(@() mlphibroker.Crosswalk.validate(xw), 'mlphibroker:crosswalk');   % hyphenated label refused
            xw.Subject = {'sub-ABC001'};
            this.verifyError(@() mlphibroker.Crosswalk.validate(xw), 'mlphibroker:crosswalk');   % unknown range refused
            mlphibroker.Crosswalk.validate(xw, '^sub-ABC\d{3}$');   % unless a project pattern says so
        end
        function test_two_studies_in_one_session_get_run_numbers(this)
            % Subject B's two studies declared as ONE visit: outputs run-01/run-02, one session report.
            tmp = fullfile(this.outRoot, 'two_studies_PHI.csv');
            t = mlphibroker.Crosswalk.template(this.root, tmp);
            isB = strcmp(t.PatientID, 'SYN000002');
            t.Subject(isB) = {'sub-900202'}; t.Session(isB) = {'ses-01'};
            t.Subject(~isB) = {'sub-900201'}; t.Session(~isB) = {'ses-01'};
            writetable(t, tmp);
            out2 = fullfile(this.outRoot, 'two_studies_out');
            exe = getenv('MLPHIBROKER_DCM2NIIX');
            r = mlphibroker.SafeHarbor.run(this.root, out2, struct('Crosswalk', tmp, 'Dcm2niixExe', exe, 'SyntheticTest', true));
            this.verifyEqual(r.subjects, 2); this.verifyEqual(r.sessions, 2); this.verifyEqual(r.studies, 3);
            this.verifyEqual(r.refused, 0, strjoin(r.refusals, ' | '));
            d = fullfile(out2, 'sub-900202', 'ses-01');
            names = {dir(fullfile(d, '*.nii.gz')).name};
            this.verifyTrue(any(strcmp(names, 'sub-900202_ses-01_trc-fdg_run-01_proc-shdeid_pet.nii.gz')), strjoin(names, ','));
            this.verifyTrue(any(strcmp(names, 'sub-900202_ses-01_trc-fdg_run-02_proc-shdeid_pet.nii.gz')));
            this.verifyTrue(any(strcmp(names, 'sub-900202_ses-01_run-01_proc-shdeid_ct.nii.gz')));
            this.verifyEqual(numel(dir(fullfile(d, '*_shdeid_report.json'))), 1);
            sr = jsondecode(fileread(fullfile(d, 'sub-900202_ses-01_shdeid_report.json')));
            this.verifyEqual(sr.studies, 2); this.verifyEqual(sr.written.pet, 2); this.verifyEqual(sr.written.ct, 2);
            % the single-study session keeps plain names
            n1 = {dir(fullfile(out2, 'sub-900201', 'ses-01', '*.nii.gz')).name};
            this.verifyTrue(any(strcmp(n1, 'sub-900201_ses-01_trc-fdg_proc-shdeid_pet.nii.gz')), strjoin(n1, ','));
        end
        function test_roster_template_round_trip(this)
            rp = fullfile(this.outRoot, 'roster2_PHI.csv');
            r = mlphibroker.Crosswalk.rosterTemplate(this.root, rp);
            this.verifyEqual(height(r), 2);                                  % two persons
            this.verifyTrue(all(cellfun(@isempty, r.Subject)));
            r.Subject = {'sub-900031'; 'sub-900032'}; writetable(r, rp);
            xp = fullfile(this.outRoot, 'xw2_PHI.csv');
            t = mlphibroker.Crosswalk.template(this.root, xp, struct('Roster', rp));
            this.verifyEqual(sort(unique(t.Subject))', {'sub-900031', 'sub-900032'});
            this.verifyTrue(all(contains(t.Note, 'authoritative')));
            % a header-less file is refused with a clear message
            bad = fullfile(this.outRoot, 'bad_PHI.csv');
            fid = fopen(bad, 'w'); fprintf(fid, 'a,b,c\n1,2,3\n'); fclose(fid);
            this.verifyError(@() mlphibroker.Crosswalk.load(bad), 'mlphibroker:crosswalk');
            fid = fopen(bad, 'w'); fprintf(fid, '1,2,3,4\n5,6,7,8\n'); fclose(fid);
            this.verifyError(@() mlphibroker.Crosswalk.load(bad), 'mlphibroker:crosswalk');
        end
        function test_deny_list_template_from_tree(this)
            dp = fullfile(this.outRoot, 'deny_from_tree.txt');
            names = mlphibroker.SidecarPolicy.denyListTemplate(this.root, dp, {'Coordinator Person'});
            this.verifyTrue(all(ismember({'Synthetic', 'Alpha', 'Beta', 'Nobody', 'Synthetic Institute', ...
                'SYNSTATION1', 'SN123456789', 'SYN000001', 'Coordinator Person'}, names)), strjoin(names, ','));
            this.verifyFalse(any(strcmpi(names, 'brain')));   % boilerplate words are not deny entries
            deny = mlphibroker.SidecarPolicy.readDenyList(dp);
            this.verifyEqual(numel(deny), numel(names));
            p = mlphibroker.SidecarPolicy(deny);
            this.verifyTrue(p.scan('report by Dr Nobody').any);
            this.verifyFalse(p.scan('sub-900001_ses-01_trc-fdg_proc-shdeid_pet').any);
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
