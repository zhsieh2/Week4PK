% Parameters for the two-compartment model for Moxifloxacin 
%https://www.thepharmajournal.com/archives/2022/vol11issue6S/PartAD/S-11-6-374-761.pdf
Vc = 0.42;   % Volume of distribution in central compartment (liters/kg)
VB = 1.31;   % Volume of distribution in peripheral compartment (liters/kg)
k12 = 3.41;  % Rate constant for central --> peripheral (h^-1)
k21 = 1.80;  % Rate constant for peripheral --> central (h^-1)
ke = 0.25;   % Rate constant for elimination (h^-1)

% Simulation parameters
tspan = [0 24];  % Simulation time (hours)
initial_conditions = [300 0];  
% Initial drug mass in central and peripheral compartments (mg)

% Solve the differential equations using ode45
[t, drug_amounts] = ode45(@(t, y) ode_equations(y, Vc, VB, k12, k21, ke), tspan, initial_conditions);

% Plot drug amounts in the central and peripheral compartments
figure;
plot(t, drug_amounts(:, 1), 'r-', t, drug_amounts(:, 2), 'b-');
xlabel('Time (hours)');
ylabel('Drug Amount (mg)');
legend('Central Compartment', 'Peripheral Compartment');
title('Drug Amount vs. Time (Two-Compartment Model)for Moxifloxacin ');

function dydt = ode_equations(y, Vc, VB, k12, k21, ke)
    % Differential equations for the two-compartment model
    Cc = y(1) / Vc;  % Concentration in the central compartment (mg/l)
    Cp = y(2) / VB;  % Concentration in the peripheral compartment (mg/l)

    % Rate of change of drug amounts in each compartment
    dydt = [k21 * Cp - k12 * Cc - ke * Cc;  % Rate of change in the central compartment
            k12 * Cc - k21 * Cp  % Rate of change in the peripheral compartment
            ]; 
end