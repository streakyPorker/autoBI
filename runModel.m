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

    names = metabolites(1,:);
    output = cell2table(cell(size(metabolites, 2), length(prop_vec)), 'VariableNames', prop_vec);
    output.Properties.RowNames = names;

    for i = 1:height(reactions)
        % calc dCdt
        concr_vec = find(strcmp(names,))
        J_plus = flux(table2array(reactions(i, 'kp')), CE, )

    end

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
