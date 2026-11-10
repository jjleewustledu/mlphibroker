classdef Test_TracerTable < matlab.unittest.TestCase
    %% The isotope is never inferred: name AND half-life must agree.
    %
    %  Created 09-Nov-2026 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlphibroker/test/+mlphibroker_unittest.
    %  Developed on Matlab 26.1 (R2026a). Copyright 2026 John J. Lee.

    methods (Test)
        function test_fdg(this)
            [tok, iso, hl] = mlphibroker.TracerTable.resolve(struct('Radiopharmaceutical', 'Fludeoxyglucose F18', 'RadionuclideHalfLife', 6586.2));
            this.verifyEqual(tok, 'fdg'); this.verifyEqual(iso, '18F'); this.verifyEqual(hl, 6586.2);
        end
        function test_co_synonym_and_water(this)
            [tok, iso] = mlphibroker.TracerTable.resolve(struct('Radiopharmaceutical', 'CO', 'RadionuclideHalfLife', 122.24));
            this.verifyEqual(tok, 'co'); this.verifyEqual(iso, '15O');
            tok = mlphibroker.TracerTable.resolve(struct('Radiopharmaceutical', '[15O]carbon monoxide', 'RadionuclideHalfLife', 122.24));
            this.verifyEqual(tok, 'co');
            tok = mlphibroker.TracerTable.resolve(struct('Radiopharmaceutical', 'Water O-15', 'RadionuclideHalfLife', 122.24));
            this.verifyEqual(tok, 'ho');
            tok = mlphibroker.TracerTable.resolve(struct('Radiopharmaceutical', 'Oxygen O15', 'RadionuclideHalfLife', 122.24));
            this.verifyEqual(tok, 'oo');
        end
        function test_bids_pet_tracername_key(this)
            % dcm2niix >= v1.0.20261029 writes TracerName instead of Radiopharmaceutical.
            [tok, iso] = mlphibroker.TracerTable.resolve(struct('TracerName', 'Fludeoxyglucose F18', 'RadionuclideHalfLife', 6586.2));
            this.verifyEqual(tok, 'fdg'); this.verifyEqual(iso, '18F');
        end
        function test_halflife_disagreement_refused(this)
            this.verifyError(@() mlphibroker.TracerTable.resolve(struct('Radiopharmaceutical', 'Fludeoxyglucose F18', 'RadionuclideHalfLife', 122.24)), 'mlphibroker:isotope');
        end
        function test_missing_halflife_refused(this)
            this.verifyError(@() mlphibroker.TracerTable.resolve(struct('Radiopharmaceutical', 'Fludeoxyglucose F18')), 'mlphibroker:halflife');
        end
        function test_unknown_name_refused_unless_override_agrees(this)
            sc = struct('Radiopharmaceutical', 'mystery compound', 'RadionuclideHalfLife', 6586.2);
            this.verifyError(@() mlphibroker.TracerTable.resolve(sc), 'mlphibroker:tracer');
            tok = mlphibroker.TracerTable.resolve(sc, struct('Tracer', 'fdg')); this.verifyEqual(tok, 'fdg');
            this.verifyError(@() mlphibroker.TracerTable.resolve(sc, struct('Tracer', 'ho')), 'mlphibroker:isotope');
        end
    end
end
