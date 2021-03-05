function [t, CR,eqns] = runModel(data_path, TC, CEx, C0, tspan)
    % Inputs:
    %   T: Temperature (Celcius)
    %   CEx: Extracellular concentrations of all metabolites (0 for some)
    %   C0: Initial concentrations of all intracellular metabolites
    %   tspan: Time to run the model. If to a steady state, tspan = -1
    % Outputs:
    %   t: Time. If tspan = -1, t = -1
    %   C: Concentration matrix.
    %   Jf: Forward fluxes associated with each of the reactions (account for
    %   reaction stoichiometry)
    %   Jb: Back fluxes associated with each of the reactions (account for
    %   reaction stoichiometry)
    %   anet: Net isotopic fractionation associated with each of the reaction



    % data_path = 'output/setup.mat';
    load(data_path,'metabolites', 'transport', 'reactions');

    R_sym = sym('R', [size(metabolites, 2), 1]);
    C_sym = sym('C', [size(metabolites, 2), 1]);

    TK = TC + 273.15;
    % reactions{i + 1, name_map('km')} = reactions{i + 1, name_map('kp')} / Keq(reactions{i + 1, name_map('DGro')}
    for i = 1:height(reactions)
        reactions.km(i) = reactions.kp(i) / Keq(reactions.DGro(i), TK);
    end

    prop_vec = {'dCdt', 'dCRdt', 'dRdt'};
    CE = 1; % need to change
    R0 = cell2mat(metabolites(2, :));

    names = metabolites(1, :);
    eqns = sym('tmp', [size(metabolites, 2), length(prop_vec)]);

    for i = 1:height(reactions)

        r_vec = reactions.r_vec{i};
        p_vec = reactions.p_vec{i};
        r_index = [];
        p_index = [];

        for j = 1:size(r_vec, 1)
            r_index(end + 1) = find(strcmp(names, r_vec(j)));
        end

        for j = 1:size(p_vec, 1)
            p_index(end + 1) = find(strcmp(names, p_vec(j)));
        end

        % calculate the flux of both direction
        J_plus = flux(reactions.kp(i), CE, C_sym(r_index), reactions.Kr_vec{i}, reactions.nr_vec{i}, ...
            C_sym(p_index), reactions.Kp_vec{i}, reactions.np_vec{i});

        J_minus = flux(reactions.km(i), CE, C_sym(p_index), reactions.Kp_vec{i}, reactions.np_vec{i}, ...
            C_sym(r_index), reactions.Kr_vec{i}, reactions.nr_vec{i});

        % calculate dCdt & dCRdt
        for j = r_index
            eqns(j, 1) = eqns(j, 1) + J_plus - J_minus;

            for k = p_index
                eqns(j, 2) = eqns(j, 2) - J_plus * R_sym(j) * reactions.ap(i) + J_minus * R_sym(k) * reactions.am(i);
            end

        end

        for j = p_index
            eqns(j, 1) = eqns(j, 1) - J_plus + J_minus;

            for k = r_index
                eqns(j, 2) = eqns(j, 2) + J_plus * R_sym(k) * reactions.ap(i) - J_minus * R_sym(j) * reactions.am(i);
            end

        end
    end

    % calc the transport position
    for i = 1:size(transport, 1)
        pos = find(strcmp(names, transport{i, 1}));
        dict = transport{i, 2};

        if dict("Mode") == "Diffusion"
            P = dict("P");
            J = -P * (C_sym(pos) - CEx(pos));
            in = (J / abs(J) + 1) / 2;
            out = (J / abs(J) - 1) / 2;
            
            eqns(pos, 1) = eqns(pos, 1) + J;
            eqns(pos, 2) = eqns(pos, 2) + J * R0(pos) * in + J * R_sym(pos) * out;
        end
    end

    % calculate dRdt
    for i = 1:size(eqns, 1)
        eqns(i, 3) = (eqns(i, 2) - R_sym(i) * eqns(i, 1)) / C_sym(i);
    end

    eqns = subs(eqns, sym('tmp', size(eqns)), zeros(size(eqns)));

    options = odeset('NonNegative', 1:length(C0) * 2);

    R0 = cell2mat(metabolites(2, :)');
    ini_val = [C0; R0];

    if tspan > 0

        [t, CR] = ode15s(@(t, CR) all_eqns(t, CR, eqns, C_sym, R_sym), [0, tspan], ini_val, options);
    else
        % Run the model until a steady state is achieved
        atSS = 0;
        tspan = 3600;

        while atSS == 0

            % Run the model with a starting tspan
            [t, CR] = ode15s(@(t, CR) all_eqns(t, CR, eqns, C_sym, R_sym), [0, tspan], ini_val, options);
            t(end, :)
            CR(end - 30:end, :)
            % Check for convergence to a steady state (change in C over
            %   2nd half of simulation is less than a specified tolerance)
            convTol = 1e-4;
            ihalf = find(t <= t(end) / 2);
            ihalf = ihalf(end);
            % ihalf = find(abs(t - t(end)) == min(abs(t - t(end)))); ihalf = ihalf(1);

            % If no convergence, increase tspan
            if any(abs(CR(ihalf, :) - CR(end, :)) ./ CR(end, :) > convTol)
                abs(CR(ihalf, :) - CR(end, :)) ./ CR(end, :)
                tspan = tspan * 2;
            else
                atSS = 1;
            end

        end

    end

    % {t,C}
    % plot(t, C);

end

function eqn = all_eqns(t, CR, eqns, C_sym, R_sym)
    eqn = [eqns(:, 1); eqns(:, 1)];
    eqn = subs(eqn, C_sym, CR(1:length(C_sym)));
    eqn = subs(eqn, R_sym, CR(length(R_sym) + 1:2 * length(R_sym)));
    eqn = eval(eqn);

end

% function eqns = processInfo(metabolites, transport, reactions, CEx, R_sym, C_sym)
% end

function output = Keq(DGro, T)
    R = 8.314; %J/(mol*K)
    output = exp(-DGro / (R * T));
end

function J = flux(k, E, concr_vec, Kr_vec, nr_vec, concp_vec, Kp_vec, np_vec)
    coef = k * E;
    numerator = prod((concr_vec ./ Kr_vec).^nr_vec);
    denominator = 1 + numerator + prod((concp_vec ./ Kp_vec).^np_vec);
    J = coef * numerator / denominator;
end
