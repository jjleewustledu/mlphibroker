classdef TracerTable
    %TRACERTABLE  trc-<token> from the radiopharmaceutical NAME and the DICOM half-life, jointly.
    %   [token, isotope, halflife_s] = mlphibroker.TracerTable.resolve(sidecar, opts)
    %   The isotope is never inferred (SPEC-PHIBROKER §4.7 item 4; idif2026 KI-66): the name must
    %   match a table row AND RadionuclideHalfLife must agree with that row's isotope within 2 %,
    %   or the series is refused. opts.Tracer may name the token when the name is ambiguous; it
    %   is still checked against the half-life. Tokens follow config/phibroker/tracers.yaml;
    %   half-lives (s) follow idif2026.PETUtilities (Turku decay table).

    properties (Constant)
        ROWS = { ...  token, isotope, halflife_s (idif2026.PETUtilities values), name patterns (regex)
            'ho',   '15O',  122.2416,     {'water', '\bh2[\s-]?15?o\b', '\[15o\]h2o'}; ...
            'oo',   '15O',  122.2416,     {'oxygen(?!.*(monoxide|dioxide))', '\bo2\b', '\[15o\]o2'}; ...
            'co',   '15O',  122.2416,     {'carbon monoxide', '(?<![a-z])co(?![a-z0-9])'}; ...
            'fdg',  '18F',  109.771*60,   {'fludeoxyglucose', 'fluorodeoxyglucose', '\bfdg\b'}; ...
            'pib',  '11C',  20.340253*60, {'\bpib\b', 'pittsburgh'}; ...
            'raclopride', '11C', 20.340253*60, {'raclopride'}; ...
            % DOTA-ECL1i is labelled with 68Ga or 64Cu (operator ruling 2026-11-09, idif2026 KI-66):
            % the same name patterns, told apart ONLY by the DICOM half-life.
            'gadotaecl1i', '68Ga', 67.719*60,  {'ecl1i', 'dota'}; ...
            'cudotaecl1i', '64Cu', 12.7003*3600, {'ecl1i', 'dota'}; ...
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
            if isempty(rows)
                error('mlphibroker:tracer', ['radiopharmaceutical name matched 0 tracer rows; ' ...
                    'pass opts.Tracer explicitly (isotope is never inferred)']);
            end
            % Several rows share a name (DOTA-ECL1i): the DICOM half-life decides among them.
            agree = rows(abs(halflife_s - cell2mat(T(rows, 3))) <= 0.02 * cell2mat(T(rows, 3)));
            if isempty(agree)
                error('mlphibroker:isotope', ['RadionuclideHalfLife %.1f s disagrees with every tracer row ' ...
                    'the name matched (%s): the name and the half-life must agree'], halflife_s, strjoin(T(rows, 1)', ', '));
            end
            if numel(agree) > 1
                error('mlphibroker:tracer', ['radiopharmaceutical name and half-life leave %d tracer rows (%s); ' ...
                    'pass opts.Tracer explicitly'], numel(agree), strjoin(T(agree, 1)', ', '));
            end
            token = T{agree, 1}; isotope = T{agree, 2};
        end
    end
end
