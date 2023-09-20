% Parameters for the two-compartment model for Acetaminophen
%https://pubmed.ncbi.nlm.nih.gov/17202677/
V1 = 30.9;  % Distribution Volume in central compartment (l)
V2 = 30.9;  % Distribution Volume in peripheral compartment (l)
k12 = 2.4;   % Rate constant for central --> peripheral (h^-1)
k21 = 1.3; % Rate constant for peripheral --> central (h^-1)
ke = 1.3; % Rate constant for elimination (h^-1)

% Simulation parameters
tspan = [0 24];  % Simulation time (hours)
initial_conditions = [300 0];  
% Initial drug mass in central and peripheral compartments (mg)

% Solve the differential equations using ode45
[t, drug_amounts] = ode45(@(t, y) ode_equations(y, V1, V2, k12, k21, ke), tspan, initial_conditions);

% Plot drug amounts in the central and peripheral compartments
figure;
plot(t, drug_amounts(:, 1), 'r-', t, drug_amounts(:, 2), 'b-');
xlabel('Time (hours)');
ylabel('Drug Amount (mg)');
legend('Central Compartment', 'Peripheral Compartment');
title('Drug Amount vs. Time (Two-Compartment Model)for Acetaminophen');

% Set y-axis ticks to be every 50 units from 0 to 300
yticks(0:50:300);

function dydt = ode_equations(y, V1, V2, k12, k21, ke)
    % Differential equations for the two-compartment model
    C1 = y(1) / V1;  % Concentration in the central compartment (μg/ml)
    C2 = y(2) / V2;  % Concentration in the peripheral compartment (μg/ml)

    % Rate of change of drug amounts in each compartment
    dydt = [k21 * C2 - k12 * C1 - ke * C1;  % Rate of change in the central compartment
            k12 * C1 - k21 * C2  % Rate of change in the peripheral compartment
            ]; 
end

