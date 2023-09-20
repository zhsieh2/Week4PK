% Parameters for the two-compartment model for Metronidazole 
%https://pubmed.ncbi.nlm.nih.gov/7459894/
Vc = 0.41;   % Volume of distribution in central compartment (liters/kg)
VB = 1.02;   % Volume of distribution in peripheral compartment (liters/kg)
k12 = 1.18;  % Rate constant for central --> peripheral (h^-1)
k21 = 0.86;  % Rate constant for peripheral --> central (h^-1)
ke = 0.22;   % Rate constant for elimination (h^-1)

% Simulation parameters
tspan = [0 24];  % Simulation time (hours)
initial_conditions = [500 0];  
% Initial drug mass in central and peripheral compartments (mg)

% Solve the differential equations using ode45
[t, drug_amounts] = ode45(@(t, y) ode_equations(y, Vc, VB, k12, k21, ke), tspan, initial_conditions);

% Plot drug amounts in the central and peripheral compartments
figure;
plot(t, drug_amounts(:, 1), 'r-', t, drug_amounts(:, 2), 'b-');
xlabel('Time (hours)');
ylabel('Drug Amount (mg)');
legend('Central Compartment', 'Peripheral Compartment');
title('Drug Amount vs. Time (Two-Compartment Model)for Metronidazole ');

function dydt = ode_equations(y, Vc, VB, k12, k21, ke)
    % Differential equations for the two-compartment model
    Cc = y(1) / Vc;  % Concentration in the central compartment (mg/l)
    Cp = y(2) / VB;  % Concentration in the peripheral compartment (mg/l)

    % Rate of change of drug amounts in each compartment
    dydt = [k21 * Cp - k12 * Cc - ke * Cc;  % Rate of change in the central compartment
            k12 * Cc - k21 * Cp  % Rate of change in the peripheral compartment
            ]; 
end
