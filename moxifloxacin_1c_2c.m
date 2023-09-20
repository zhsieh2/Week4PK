% Define the differential equations function for the two-compartment model
function dydt = ode_equations(t, y)
    % Parameters for the two-compartment model for Moxifloxacin
    Vc = 0.42;   % Volume of distribution in central compartment (liters/kg)
    VB = 1.31;   % Volume of distribution in peripheral compartment (liters/kg)
    k12 = 3.41;  % Rate constant for central --> peripheral (h^-1)
    k21 = 1.80;  % Rate constant for peripheral --> central (h^-1)
    ke = 0.25;   % Rate constant for elimination (h^-1)

    % Concentrations in the central and peripheral compartments
    Cc = y(1) / Vc;  % Concentration in the central compartment (mg/l)
    Cp = y(2) / VB;  % Concentration in the peripheral compartment (mg/l)

    % Rate of change of drug amounts in each compartment
    dydt = [
        k21 * Cp - k12 * Cc - ke * Cc; % Rate of change in the central compartment
        k12 * Cc - k21 * Cp           % Rate of change in the peripheral compartment
    ];
end

% Parameters for the one-compartment model (Moxifloxacin)
k_moxifloxacin = 0.25; % Elimination rate constant 
% Define the time span for the simulation
tspan_moxifloxacin = [0 24]; % Start at time 0 and end at time 24
% Initial drug amount
C0_moxifloxacin = 300; % Initial amount (mg)
% Solve the differential equation using ode45
[t_moxifloxacin, C_moxifloxacin] = ode45(@(t, C) (-k_moxifloxacin)*C,  tspan_moxifloxacin, C0_moxifloxacin);

% Parameters for the two-compartment model (Moxifloxacin)
Vc_moxifloxacin = 0.42;   % Volume of distribution in central compartment (liters/kg)
VB_moxifloxacin = 1.31;   % Volume of distribution in peripheral compartment (liters/kg)
k12_moxifloxacin = 3.41;  % Rate constant for central --> peripheral (h^-1)
k21_moxifloxacin = 1.80;  % Rate constant for peripheral --> central (h^-1)
ke_moxifloxacin = 0.25;   % Rate constant for elimination (h^-1)

% Define the time span for the simulation (Moxifloxacin)
tspan_moxifloxacin = [0 24]; % Start at time 0 and end at time 24
% Initial drug amount (Moxifloxacin)
C0_moxifloxacin_2c = [300 0]; % Initial amounts in central and peripheral compartments (mg)

% Solve the differential equations using ode45 (Moxifloxacin)
[t_moxifloxacin_2c, drug_amounts_moxifloxacin_2c] = ode45(@ode_equations, tspan_moxifloxacin, C0_moxifloxacin_2c);

% Plot both models on the same figure
figure;
subplot(2,1,1);
plot(t_moxifloxacin, C_moxifloxacin, 'b-');
xlabel('Time (hours)');
ylabel('Drug Amount (mg)');
title('One-Compartment Model (Moxifloxacin)');
grid on;

subplot(2,1,2);
plot(t_moxifloxacin_2c, drug_amounts_moxifloxacin_2c(:, 1), 'r-', t_moxifloxacin_2c, drug_amounts_moxifloxacin_2c(:, 2), 'g-');
xlabel('Time (hours)');
ylabel('Drug Amount (mg)');
legend('Central Compartment', 'Peripheral Compartment');
title('Two-Compartment Model (Moxifloxacin)');
grid on;
