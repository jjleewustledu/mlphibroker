classdef SafeHarbor
    %SAFEHARBOR  Batch conversion of a heterogeneous DICOM tree to Safe Harbor PET/CT/ROI files.
    %
    %   report = mlphibroker.SafeHarbor.run(root, outRoot, opts)
    %
    %   Runs on the PHI side (MATTA). Traverses root with dicomCollection, finds the unique
    %   subjects and studies, maps them to canonical codes through a crosswalk, and for every
    %   study writes under outRoot/<sub>/<ses>/:
    %     <sub>_<ses>_trc-<token>_proc-shdeid_pet.nii.gz + .json   dynamic or static PET
    %     (a session holding several studies, e.g. a dynamic and a delayed scan of one injection,
    %      gets run-01, run-02 ... in acquisition order on every PET, CT and ROI file)
    %     <sub>_<ses>_proc-shdeid_ct.nii.gz + .json                 CT
    %     <sub>_<ses>_roi-<name>_proc-shdeid_mask.nii.gz + .json    one per RTSTRUCT ROI (MIM)
    %     <sub>_<ses>_shdeid_report.json                             counts only
    %   and outRoot/shdeid_batch_report.json (counts only). Voxel data are converted by dcm2niix
    %   (PET, CT) or dicomContours/createMask (RTSTRUCT); everything identifying is removed by
    %   mlphibroker.SidecarPolicy (allow-list), mlphibroker.Nifti1 (header scrub) and the
    %   PHI-shape scan with quarantine. Times are seconds from scan start; dates never appear.
    %
    %   opts fields (all optional unless stated):
    %     Crosswalk      path to a CSV with columns PatientID,Subject,StudyInstanceUID,Session
    %                    (PHI-side file from the key vault). If absent and AutoCrosswalk is true,
    %                    synthetic codes sub-9000NN and random ses-NN are minted and written to
    %                    outRoot/crosswalk_PHI.csv, which MUST stay in the PHI zone.
    %     AutoCrosswalk  logical, default false
    %     Dcm2niixExe    path to dcm2niix(.exe); default: found on PATH
    %     DenyListFile   one staff/participant name per line (key vault)
    %     Tracer         token override when the radiopharmaceutical name is ambiguous
    %     SyntheticTest  logical; permits running under Claude Code for sub-9xxxxx fixtures only
    %
    %   Governing text: idif2026 specs/SPEC_phi_broker.md §0 item 5, §4.2, §5, §13; WO-0002.
    %   This is a fast working example for time-sensitive kinetics, not the enclave ETL.
    %
    %  Created 09-Nov-2026 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlphibroker/src/+mlphibroker.
    %  Developed on Matlab 26.1 (R2026a). Copyright 2026 John J. Lee.

    properties (Constant)
        TOOL = 'mlphibroker'
        VERSION = '0.1.1'
    end

    methods (Static)
        function report = run(root, outRoot, opts)
            arguments
                root {mustBeFolder}
                outRoot (1,:) char
                opts struct = struct()
            end
            import mlphibroker.*
            synthetic = isfield(opts, 'SyntheticTest') && isequal(opts.SyntheticTest, true);
            if ~isempty(getenv('CLAUDECODE')) && ~synthetic
                error('mlphibroker:claudecode', ...
                    'refusing to run inside a Claude Code session (SPEC-REDCAP-PHI §7)');
            end
            if ~isfolder(outRoot), mkdir(outRoot); end

            tree = DicomTree(root);
            xw = Crosswalk.resolve(tree, outRoot, opts);
            if synthetic && any(~startsWith(xw.Subject, 'sub-9'))
                error('mlphibroker:synthetic', 'SyntheticTest permits sub-9xxxxx subjects only');
            end
            deny = {};
            if isfield(opts, 'DenyListFile') && ~isempty(opts.DenyListFile)
                deny = SidecarPolicy.readDenyList(opts.DenyListFile);
            end
            policy = SidecarPolicy(deny);
            conv = Dcm2niix(opts);

            report = struct('tool', SafeHarbor.TOOL, 'version', SafeHarbor.VERSION, ...
                'dcm2niix_version', conv.version(), ...
                'subjects', numel(unique(xw.Subject)), 'sessions', height(unique(xw(:, {'Subject', 'Session'}))), 'studies', height(xw), ...
                'series_by_kind', tree.countByKind(), 'written', struct('pet', 0, 'ct', 0, 'roi', 0), ...
                'refused', 0, 'unhandled_series', 0, 'quarantined', 0, 'dropped_key_names', {{}}, ...
                'refusals', {{}});
            dropped = {};
            groups = unique(xw(:, {'Subject', 'Session'}), 'stable');
            for g = 1:height(groups)
                sub = char(groups.Subject(g)); ses = char(groups.Session(g));
                studyUids = xw.StudyInstanceUID(strcmp(xw.Subject, sub) & strcmp(xw.Session, ses));
                studyUids = tree.orderStudiesByAcquisition(studyUids);   % PHI-side; only the ORDER is kept
                outDir = fullfile(outRoot, sub, ses);
                if ~isfolder(outDir), mkdir(outDir); end
                sr = struct('subject', sub, 'session', ses, 'studies', numel(studyUids), ...
                    'written', struct('pet', 0, 'ct', 0, 'roi', 0), 'refused', 0, 'unhandled', 0, ...
                    'quarantined', 0, 'refusals', {{}});
                written = containers.Map('KeyType', 'char', 'ValueType', 'any');  % SeriesInstanceUID -> converted .nii.gz
                for st = 1:numel(studyUids)
                    rows = tree.seriesOfStudy(studyUids{st});
                    run = 0; if numel(studyUids) > 1, run = st; end   % run-NN only when a session holds several studies
                    for r = 1:height(rows)
                        kind = char(rows.Kind(r));
                        try
                            switch kind
                                case {'pet_dynamic', 'pet_static'}
                                    [n, d, q, paths] = PetSession.convert(rows(r, :), sub, ses, outDir, conv, policy, opts, run);
                                    sr.written.pet = sr.written.pet + n; dropped = union(dropped, d);
                                    sr.quarantined = sr.quarantined + q;
                                    if n > 0, written(rows.SeriesInstanceUID{r}) = paths{1}; end
                                case 'ct'
                                    [n, d, q, paths] = CtSeries.convert(rows(r, :), sub, ses, outDir, conv, policy, run);
                                    sr.written.ct = sr.written.ct + n; dropped = union(dropped, d);
                                    sr.quarantined = sr.quarantined + q;
                                    if n > 0, written(rows.SeriesInstanceUID{r}) = paths{1}; end
                                case 'rtstruct'
                                    % handled after the images so the reference geometry exists
                                otherwise
                                    sr.unhandled = sr.unhandled + 1;
                            end
                        catch err
                            sr.refused = sr.refused + 1;
                            sr.refusals{end+1} = sprintf('%s series: %s', kind, strtok(err.message, newline));
                        end
                    end
                    rs = rows(strcmp(rows.Kind, 'rtstruct'), :);
                    for r = 1:height(rs)
                        try
                            [n, q] = RtStruct.convert(rs(r, :), rows, written, sub, ses, outDir, policy, run);
                            sr.written.roi = sr.written.roi + n; sr.quarantined = sr.quarantined + q;
                        catch err
                            sr.refused = sr.refused + 1;
                            sr.refusals{end+1} = sprintf('rtstruct: %s', strtok(err.message, newline));
                        end
                    end
                end
                Report.write(fullfile(outDir, sprintf('%s_%s_shdeid_report.json', sub, ses)), sr);
                report.written.pet = report.written.pet + sr.written.pet;
                report.written.ct = report.written.ct + sr.written.ct;
                report.written.roi = report.written.roi + sr.written.roi;
                report.refused = report.refused + sr.refused;
                report.unhandled_series = report.unhandled_series + sr.unhandled;
                report.quarantined = report.quarantined + sr.quarantined;
                report.refusals = [report.refusals, sr.refusals];
            end
            report.dropped_key_names = sort(dropped);
            Report.write(fullfile(outRoot, 'shdeid_batch_report.json'), report);
            fprintf('%s %s: %d subjects, %d studies; pet %d ct %d roi %d written; %d refused, %d quarantined\n', ...
                SafeHarbor.TOOL, SafeHarbor.VERSION, report.subjects, report.studies, ...
                report.written.pet, report.written.ct, report.written.roi, report.refused, report.quarantined);
        end
    end
end
