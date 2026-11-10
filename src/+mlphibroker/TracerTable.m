classdef TracerTable
    %TRACERTABLE  trc-<token> from the radiopharmaceutical NAME and the DICOM half-life, jointly.
    %   [token, isotope, halflife_s] = mlphibroker.TracerTable.resolve(sidecar, opts)
    %   The isotope is never inferred (SPEC-PHIBROKER §4.7 item 4; idif2026 KI-66): the name must
    %   match a table row AND RadionuclideHalfLife must agree with that row's isotope within 2 %,
    %   or the series is refused. opts.Tracer may name the token when the name is ambiguous; it
    %   is still checked against the half-life. Tokens follow config/phibroker/tracers.yaml;
    %   half-lives (s) follow idif2026.PETUtilities (Turku decay table).

    properties (Constant)
        ROWS = { ...  token, isotope, halflife_s, name patterns (case-insensitive regex)
            'ho',   '15O',  122.2416,  {'water', '\bh2[\s-]?15?o\b', '\[15o\]h2o'}; ...
            'oo',   '15O',  122.2416,  {'oxygen(?!.*(monoxide|dioxide))', '\bo2\b', '\[15o\]o2'}; ...
            'co',   '15O',  122.2416,  {'carbon monoxide', '(?<![a-z])co(?![a-z0-9])'}; ...
            'fdg',  '18F',  6586.26,   {'fludeoxyglucose', 'fluorodeoxyglucose', '\bfdg\b'}; ...
            'pib',  '11C',  1223.4,    {'\bpib\b', 'pittsburgh'}; ...
            'raclopride', '11C', 1223.4, {'raclopride'}; ...
            }
    end

    methods (Static)
        function [token, isotope, halflife_s] = resolve(sc, opts)
            T = mlphibroker.TracerTable.ROWS;
            name = '';   % dcm2niix <= 2025 writes Radiopharmaceutical; >= 2026 writes BIDS-PET TracerName
            if isfield(sc, 'Radiopharmaceutical'), name = lower(char(sc.Radiopharmaceutical)); end
            if isempty(name) && isfield(sc, 'TracerName'), name = lower(char(sc.TracerName)); end
            if ~isfield(sc, 'RadionuclideHalfLife') || isempty(sc.RadionuclideHalfLife)
                error('mlphibroker:halflife', 'RadionuclideHalfLife missing: the isotope cannot be verified');
            end
            halflife_s = double(sc.RadionuclideHalfLife);
            cands = false(size(T, 1), 1);
            for r = 1:size(T, 1)
                pats = T{r, 4};
                for p = 1:numel(pats)
                    if ~isempty(regexp(name, pats{p}, 'once')), cands(r) = true; end
                end
            end
            if nargin > 1 && isfield(opts, 'Tracer') && ~isempty(opts.Tracer)
                r = find(strcmp(T(:, 1), lower(char(opts.Tracer))), 1);
                if isempty(r), error('mlphibroker:tracer', 'opts.Tracer is not in the tracer table'); end
                cands = false(size(cands)); cands(r) = true;
            end
            rows = find(cands);
            if numel(rows) ~= 1
                error('mlphibroker:tracer', ['radiopharmaceutical name matched %d tracer rows; ' ...
                    'pass opts.Tracer explicitly (isotope is never inferred)'], numel(rows));
            end
            token = T{rows, 1}; isotope = T{rows, 2}; ref = T{rows, 3};
            if abs(halflife_s - ref) > 0.02 * ref
                error('mlphibroker:isotope', ['RadionuclideHalfLife disagrees with the %s row (%s): ' ...
                    'the name and the half-life must agree'], token, isotope);
            end
        end
    end
end
