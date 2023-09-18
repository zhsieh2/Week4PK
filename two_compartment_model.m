function two_compartment_model()
    % Parameters for the two-compartment model
    V1 = 100;  % Volume of distribution in central compartment (mL)
    V2 = 100;  % Volume of distribution in peripheral compartment (mL)
    k12 = 50;   % Rate constant for central --> peripheral
    k21 = 10; % Rate constant for peripheral --> central
    ke = 10; % Rate constant for elimination

    % Simulation parameters
    tspan = [0 60];  % Simulation time (minutes)
    initial_conditions = [1000 0];  % Initial drug mass in central and peripheral compartments (ng)

    % Solve the differential equations using ode45
    [t, drug_amounts] = ode45(@(t, y) ode_equations(y, V1, V2, k12, k21, ke), tspan, initial_conditions);

    % Plot drug amounts in the central and peripheral compartments
    figure;
    plot(t, drug_amounts(:, 1), 'r-', t, drug_amounts(:, 2), 'b-');
    xlabel('Time (min)');
    ylabel('Drug Amount (ng)');
    legend('Central Compartment', 'Peripheral Compartment');
    title('Drug Amount vs. Time (Two-Compartment Model)');
end

function dydt = ode_equations(y, V1, V2, k12, k21, ke)
    % Differential equations for the two-compartment model
    C1 = y(1) / V1;  % Concentration in the central compartment (ng/mL)
    C2 = y(2) / V2;  % Concentration in the peripheral compartment (ng/mL)

    % Rate of change of drug amounts in each compartment
    dydt = [k21 * C2 * (V2/V1) - k12 * C1 - ke * C1;  % Rate of change in the central compartment
            k12 * C1 * (V1/V2) - k21 * C2  % Rate of change in the peripheral compartment
            ]; 
end
