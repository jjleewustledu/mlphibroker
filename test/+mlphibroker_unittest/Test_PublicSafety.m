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
        function test_fixture_dicoms_are_fabricated(this)
            % Every DICOM in the committed fixture must carry a fabricated identity: this is what
            % makes the .gitignore negation for that folder safe in a public repository.
            root = fullfile(fileparts(mfilename('fullpath')), '..', 'fixtures', 'synthetic_dicom');
            files = dir(fullfile(root, '**', '*.dcm'));
            this.assertGreaterThan(numel(files), 50, 'the fixture must be present');
            for i = 1:numel(files)
                h = dicominfo(fullfile(files(i).folder, files(i).name), 'UseDictionaryVR', true);
                this.verifyTrue(startsWith(char(h.PatientName.FamilyName), 'Synthetic'), files(i).name);
                this.verifyTrue(startsWith(char(h.PatientID), 'SYN'), files(i).name);
                this.verifyTrue(startsWith(char(h.StudyDate), '2001'), files(i).name);
                this.verifyEqual(char(h.InstitutionName), 'Synthetic Institute', files(i).name);
            end
        end
    end
end
