classdef Test_SidecarPolicy < matlab.unittest.TestCase
    %% Allow-list keeps decay and timing keys, drops identity and absolute times; scan finds shapes.
    %
    %  Created 09-Nov-2026 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlphibroker/test/+mlphibroker_unittest.
    %  Developed on Matlab 26.1 (R2026a). Copyright 2026 John J. Lee.

    methods (Test)
        function test_allow_drops_identity_keeps_decay(this)
            p = mlphibroker.SidecarPolicy();
            src = struct('PatientName', 'x', 'AcquisitionTime', '12:00:11', 'TimeZero', '12:00:11', ...
                'StationName', 's', 'DeviceSerialNumber', 'd', 'InstitutionName', 'i', ...
                'DecayCorrection', 'START', 'RadionuclideHalfLife', 6586.2, 'FrameTimesStart', [0 30], ...
                'FrameDuration', [30 30], 'PatientWeight', 70.5, 'Units', 'BQML');
            [kept, dropped] = p.allow(src);
            this.verifyEqual(sort(fieldnames(kept))', sort({'DecayCorrection', 'RadionuclideHalfLife', ...
                'FrameTimesStart', 'FrameDuration', 'PatientWeight', 'Units'}));
            this.verifyEqual(sort(dropped), sort({'PatientName', 'AcquisitionTime', 'TimeZero', ...
                'StationName', 'DeviceSerialNumber', 'InstitutionName'}));
        end
        function test_allow_list_never_contains_identity_or_clock_keys(this)
            forbidden = {'PatientName', 'PatientID', 'PatientBirthDate', 'PatientSex', 'PatientAge', ...
                'AcquisitionTime', 'AcquisitionDateTime', 'AcquisitionDate', 'TimeZero', 'StudyDate', ...
                'RadiopharmaceuticalStartTime', 'StationName', 'DeviceSerialNumber', 'InstitutionName', ...
                'AccessionNumber', 'StudyInstanceUID', 'SeriesInstanceUID', 'PatientSize', 'ImageComments'};
            this.verifyEmpty(intersect(forbidden, mlphibroker.SidecarPolicy.ALLOW));
            required = {'DecayCorrection', 'RadionuclideHalfLife', 'Radiopharmaceutical', 'FrameTimesStart', ...
                'FrameDuration', 'FrameReferenceTime', 'Units', 'PatientWeight'};
            this.verifyEqual(sort(intersect(required, mlphibroker.SidecarPolicy.ALLOW)), sort(required));
        end
        function test_scan_shapes(this)
            p = mlphibroker.SidecarPolicy({'Zorbwick'});
            h = p.scan('clean text FrameTimesStart [0,30,60] dose 185000000');
            this.verifyFalse(h.any);
            h = p.scan('2001-02-03'); this.verifyEqual(h.date, 1);
            h = p.scan('20010203'); this.verifyEqual(h.date, 1);
            h = p.scan('2/3/2001'); this.verifyEqual(h.date, 1);
            h = p.scan('12:00:11.000000'); this.verifyEqual(h.clock_time, 1);
            dated = ['sub-900001_ses-' '2001' '0203' '120011'];   % built at run time: the source stays clean
            h = p.scan(dated); this.verifyEqual(h.dated_session, 1);
            h = p.scan('Dr ZORBWICK signed'); this.verifyEqual(h.deny_list, 1);
            h = p.scan('zorbwickian'); this.verifyEqual(h.deny_list, 0);
            h = p.scan('ses-01 trc-fdg proc-shdeid'); this.verifyFalse(h.any);
        end
        function test_report_refuses_dates(this)
            tmp = [tempname '.json'];
            this.verifyError(@() mlphibroker.Report.write(tmp, struct('note', '2001-02-03')), 'mlphibroker:report');
            mlphibroker.Report.write(tmp, struct('count', 3)); this.verifyTrue(isfile(tmp)); delete(tmp);
        end
    end
end
