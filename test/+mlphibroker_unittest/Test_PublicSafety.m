classdef Test_PublicSafety < matlab.unittest.TestCase
    %% The repository is public: no real subject code, dated session label, clock time or e-mail
    %  may appear in any tracked text file. The synthetic DICOM fixture is binary and fabricated.
    %
    %  Created 09-Nov-2026 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlphibroker/test/+mlphibroker_unittest.
    %  Developed on Matlab 26.1 (R2026a). Copyright 2026 John J. Lee.

    methods (Test)
        function test_tracked_text_is_clean(this)
            root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
            files = [dir(fullfile(root, '**', '*.m')); dir(fullfile(root, '**', '*.md')); ...
                dir(fullfile(root, '**', '*.py')); dir(fullfile(root, '**', '*.csv')); dir(fullfile(root, '**', '*.json'))];
            files = files(~contains({files.folder}, [filesep '.git']));
            this.assertGreaterThan(numel(files), 5, 'the scan must see the repository');
            bad = {};
            for i = 1:numel(files)
                txt = fileread(fullfile(files(i).folder, files(i).name));
                if ~isempty(regexp(txt, 'sub-108\d{3}(?!x)', 'once')), bad{end+1} = [files(i).name ': real subject code']; end %#ok<AGROW>
                if ~isempty(regexp(txt, 'ses-\d{8,}', 'once')), bad{end+1} = [files(i).name ': dated session label']; end %#ok<AGROW>
                if ~isempty(regexp(txt, '[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}', 'once')), bad{end+1} = [files(i).name ': e-mail']; end %#ok<AGROW>
            end
            this.verifyEmpty(bad, strjoin(bad, newline));
        end
    end
end
