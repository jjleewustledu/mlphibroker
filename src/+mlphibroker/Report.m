classdef Report
    %REPORT  Counts-only JSON reports (SPEC-PHIBROKER §13.3: no value, date, time, name, id, path).
    methods (Static)
        function write(path, s)
            txt = jsonencode(s, 'PrettyPrint', true);
            p = mlphibroker.SidecarPolicy();
            hits = p.scan(txt);
            if hits.date + hits.clock_time + hits.dated_session > 0
                error('mlphibroker:report', 'a report would carry a date or clock-time shape; refusing to write it');
            end
            fid = fopen(path, 'w'); fwrite(fid, txt, 'char'); fclose(fid);
        end
    end
end
