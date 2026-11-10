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
        function test_dota_ecl1i_resolved_by_halflife_only(this)
            % Same name, two nuclides: the DICOM half-life decides; a third half-life is refused.
            [tok, iso] = mlphibroker.TracerTable.resolve(struct('TracerName', 'DOTA-ECL1i', 'RadionuclideHalfLife', 4063.1));
            this.verifyEqual(tok, 'gadotaecl1i'); this.verifyEqual(iso, '68Ga');
            [tok, iso] = mlphibroker.TracerTable.resolve(struct('Radiopharmaceutical', '68Ga DOTA ECL1i', 'RadionuclideHalfLife', 4063.1));
            this.verifyEqual(tok, 'gadotaecl1i'); this.verifyEqual(iso, '68Ga');
            [tok, iso] = mlphibroker.TracerTable.resolve(struct('Radiopharmaceutical', 'Cu64-DOTA-ECL1i', 'RadionuclideHalfLife', 45721.1));
            this.verifyEqual(tok, 'cudotaecl1i'); this.verifyEqual(iso, '64Cu');
            this.verifyError(@() mlphibroker.TracerTable.resolve(struct('TracerName', 'DOTA-ECL1i', 'RadionuclideHalfLife', 122.24)), 'mlphibroker:isotope');
            % an explicit but WRONG override is refused by the half-life, never trusted
            this.verifyError(@() mlphibroker.TracerTable.resolve(struct('TracerName', 'DOTA-ECL1i', 'RadionuclideHalfLife', 4063.1), struct('Tracer', 'cudotaecl1i')), 'mlphibroker:isotope');
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
