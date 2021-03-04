function [t, C, Jf, Jb, anet] = runModel(data_path, TC, CEx, C0, tspan)
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

    data_path = 'output/setup.mat';
    load(data_path);
    R_sym = sym('R', [size(metabolites, 2), 1]);
    C_sym = sym('C', [size(metabolites, 2), 1]);

    TK = TC + 273.15;

    processInfo(metabolites, transport, reactions, R, CEx, C0, Rsym, C_sym, TK)

    options = odeset('NonNegative', 1:length(C0));

    if tspan > 0
        [t, C] = ode15s(@(t, C) all_eqns(t, C, CEx, TK), [0, tspan], C0, options);
    else
        % Run the model until a steady state is achieved
        atSS = 0;
        tspan = 3600;

        while atSS == 0
            % Run the model with a starting tspan
            [t, C] = ode15s(@(t, C) all_eqns(t, C, CEx, TK), [tspan, 0], C0, options);

            % Check for convergence to a steady state (change in C over
            %   2nd half of simulation is less than a specified tolerance)
            convTol = 1e-4;
            ihalf = find(abs(t - t(end)) == min(abs(t - t(end)))); ihalf = ihalf(1);

            % If no convergence, increase tspan
            if any(C(ihalf, :) ./ C(end, :) > convTol)
                tspan = tspan * 2;
            else
                atSS = 1;
            end

        end

    end

end

function output = processInfo(metabolites, transport, reactions, R, CEx, C0, Rsym, C_sym, TK)
    prop_vec = {'dCdt', 'dCRdt', 'dRdt'};
    CE = 1; % need to change

    names = metabolites(1, :);
    output = sym('temp',[size(metabolites, 2), length(prop_vec)]);

    for i = 1:height(reactions)

        r_vec = reactions.r_vec{i};
        p_vec = reactions.p_vec{i};
        % [r_index, p_index] = findPos(names, r_vec, p_vec);
        r_index = [];
        p_index = [];

        for j = 1:size(r_vec, 1)
            r_index(end + 1) = find(strcmp(names, r_vec(j)));
        end

        for j = 1:size(p_vec, 1)
            p_index(end + 1) = find(strcmp(names, p_vec(j)));
        end

        J_plus = flux(reactions.kp(i), CE, C_sym(r_index), reactions.Kr_vec{i}, reactions.nr_vec{i}, ...
            C_sym(p_index), reactions.Kp_vec{i}, reactions.np_vec{i});

        J_minus = flux(reactions.km(i), CE, C_sym(p_index), reactions.Kp_vec{i}, reactions.np_vec{i}, ...
            C_sym(r_index), reactions.Kr_vec{i}, reactions.nr_vec{i});

        % calc dCdt & dCRdt
        for j = r_index
            output(j, 1) = output(j, 1) + J_plus - J_minus;

            for k = p_index
                output(j,2) = output(j,2) - J_plus * R_sym(j) * reactions.ap(i) + J_minus * R_sym(k) * reactions.am(i);
            end

        end

        for j = p_index
            output(j, 1)  = output(j, 1)  - J_plus + J_minus;

            for k = r_index
                output(j,2) = output(j,2) + J_plus * R_sym(k) * reactions.ap(i) - J_minus*R_sym(j) * reactions.am(i);
            end

        end

        % calc dCRdt
    end
    for i=1:size(output,1)
        output(i,3)=(output(i,2)-R_sym(i)*output(i,1))/C_sym(i);
    end
    output=subs(output, sym('temp',size(output)), zeros(size(output)));
end





    function dCdt = dCdt_eqns(t, C, J)
        % C: Vector of metabolites (A,B,C,D)
        dCdt = zeros(size(C));

        % Specify the differential equations
        dCdt(1) = JAtrans - JA_B + JB_A;
        dCdt(2) = JA_B - JB_A - JB_CD + JCD_B;

    end

    function dRdt = dRdt_eqns(t, R, J)
        dRdt = zeros(size(R));

    end

    function J = flux(k, E, concr_vec, Kr_vec, nr_vec, concp_vec, Kp_vec, np_vec)
        coef = k * E;
        numerator = prod((concr_vec ./ Kr_vec).^nr_vec);
        denominator = 1 + numerator + prod((concp_vec ./ Kp_vec).^np_vec);
        J = coef * numerator / denominator;
    end
