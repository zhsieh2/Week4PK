function two_compartment_model()
    % Parameters for the two-compartment model
    V1 = 100;  % Volume of distribution in central compartment (mL)
    V2 = 200;  % Volume of distribution in peripheral compartment (mL)
    CL = 10;   % Clearance (mL/min)

    % Simulation parameters
    tspan = [0 60];  % Simulation time (minutes)
    initial_conditions = [1000 0 0];  % Initial drug amounts in central, peripheral, and eliminated compartments (ng)

    % Solve the differential equations using ode45
    [t, drug_amounts] = ode45(@(t, y) ode_equations(y, V1, V2, CL), tspan, initial_conditions);
    % y is the concentration, CL is the clearance

    % Plot drug amounts in the central and peripheral compartments
    figure;
    plot(t, drug_amounts(:, 1), 'r-', t, drug_amounts(:, 2), 'b-');
    xlabel('Time (min)');
    ylabel('Drug Amount (ng)');
    legend('Central Compartment', 'Peripheral Compartment');
    title('Drug Amount vs. Time (Two-Compartment Model)');
end

function dydt = ode_equations(y, V1, V2, CL)
    % Differential equations for the two-compartment model
    C1 = y(1) / V1;  % Concentration in the central compartment (ng/mL)
    C2 = y(2) / V2;  % Concentration in the peripheral compartment (ng/mL)

    % Rate of change of drug amounts in each compartment
    dydt = [-CL * C1;  % Rate of change in the central compartment
            CL * C1 - CL * C2;  % Rate of change in the peripheral compartment
            CL * C2];  % Rate of change in the eliminated compartment
end
