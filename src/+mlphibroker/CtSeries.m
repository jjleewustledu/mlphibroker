classdef CtSeries
    %CTSERIES  One CT series -> Safe Harbor .nii.gz + .json (dcm2niix, allow-list, header scrub).
    %   [n, dropped, quarantined] = mlphibroker.CtSeries.convert(row, sub, ses, outDir, conv, policy)

    methods (Static)
        function [n, dropped, quarantined, paths] = convert(row, sub, ses, outDir, conv, policy, run)
            import mlphibroker.*
            if nargin < 7 || isempty(run), run = 0; end
            n = 0; quarantined = 0; dropped = {}; paths = {};
            prefix = sprintf('%s_%s_ser-', sub, ses);
            tmpOut = tempname(); mkdir(tmpOut);
            cleanup = onCleanup(@() rmdir(tmpOut, 's'));
            [jsons, niis] = conv.convertSeries(row.Files{1}, prefix, tmpOut);
            for k = 1:numel(jsons)
                src = jsondecode(fileread(jsons{k}));
                [kept, d] = policy.allow(src); dropped = union(dropped, d);
                sc = struct('Subject', sub, 'Session', ses, 'Modality', 'CT');
                f = fieldnames(kept);
                for i = 1:numel(f), if ~strcmp(f{i}, 'Modality'), sc.(f{i}) = kept.(f{i}); end, end
                sc.SafeHarbor = struct('tool', SafeHarbor.TOOL, 'version', SafeHarbor.VERSION, ...
                    'work_order', 'WO-0002r0', 'policy', 'allow_list', 'calendar_dates', 'none', ...
                    'clock_times', 'none', 'nifti_header_scrubbed', true);
                base = sprintf('%s_%s%s_proc-shdeid_ct', sub, ses, PetSession.runTag(run, k, numel(jsons)));
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
    end
end
