classdef Dcm2niix
    %DCM2NIIX  Run dcm2niix on ONE series with the WO-0002 flags; return the .json/.nii.gz pairs.
    %   conv = mlphibroker.Dcm2niix(opts); [jsons, niis] = conv.convertSeries(files, prefix, outDir)
    %   Flags: -b y (sidecar) -ba y (anonymised BIDS) -z y (gzip) -t n (no text notes) -v 0, and
    %   -f "<prefix>%s" so a file name carries only the canonical codes and the series number
    %   (never %t date-time, %n name, %i patient id). The series' files are copied to a private
    %   temporary folder first, so a mixed export folder converts one series at a time.
    %   dcm2niix is a converter, not a de-identifier: SidecarPolicy and Nifti1 finish the job.

    properties (SetAccess = private)
        Exe
    end

    methods
        function this = Dcm2niix(opts)
            if nargin > 0 && isfield(opts, 'Dcm2niixExe') && ~isempty(opts.Dcm2niixExe)
                this.Exe = char(opts.Dcm2niixExe);
            else
                this.Exe = mlphibroker.Dcm2niix.findOnPath();
            end
            if isempty(this.Exe) || ~isfile(this.Exe)
                error('mlphibroker:dcm2niix', 'dcm2niix executable not found; pass opts.Dcm2niixExe');
            end
        end

        function [jsons, niis, log] = convertSeries(this, files, prefix, outDir)
            files = cellstr(files);
            tmpIn = tempname(); mkdir(tmpIn);
            tmpOut = tempname(); mkdir(tmpOut);
            cleanup = onCleanup(@() mlphibroker.Dcm2niix.rmdirs({tmpIn, tmpOut}));
            for i = 1:numel(files)
                copyfile(files{i}, fullfile(tmpIn, sprintf('%06d.dcm', i)));
            end
            cmd = sprintf('"%s" -b y -ba y -z y -t n -v 0 -w 1 -f "%s%%s" -o "%s" "%s"', ...
                this.Exe, prefix, tmpOut, tmpIn);
            [rc, log] = system(cmd);
            if rc ~= 0
                error('mlphibroker:dcm2niix', 'dcm2niix exit code %d', rc);
            end
            js = dir(fullfile(tmpOut, [prefix '*.json']));
            if isempty(js), error('mlphibroker:dcm2niix', 'dcm2niix produced no sidecar'); end
            if ~isfolder(outDir), mkdir(outDir); end
            jsons = {}; niis = {};
            for i = 1:numel(js)
                j = fullfile(js(i).folder, js(i).name);
                n = regexprep(j, '\.json$', '.nii.gz');
                if ~isfile(n), continue; end
                movefile(j, fullfile(outDir, js(i).name));
                movefile(n, fullfile(outDir, regexprep(js(i).name, '\.json$', '.nii.gz')));
                jsons{end+1} = fullfile(outDir, js(i).name); %#ok<AGROW>
                niis{end+1} = fullfile(outDir, regexprep(js(i).name, '\.json$', '.nii.gz')); %#ok<AGROW>
            end
        end
    end

    methods (Static)
        function exe = findOnPath()
            exe = '';
            if ispc, [rc, out] = system('where dcm2niix.exe'); else, [rc, out] = system('which dcm2niix'); end
            if rc == 0
                lines = strtrim(strsplit(strtrim(out), newline));
                exe = lines{1};
            end
        end

        function rmdirs(dirs)
            for i = 1:numel(dirs)
                if isfolder(dirs{i}), rmdir(dirs{i}, 's'); end
            end
        end
    end
end
