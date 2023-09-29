% Parameters for the one-compartment model with plasma protein binding
k_e = 0.1;          % Elimination rate constant
k_association = 30;  % Association rate constant for plasma protein binding
k_dissociation = 20; % Dissociation rate constant for plasma protein binding
V_one_compartment = 1;  % Volume of distribution for one-compartment model
F = 0.8;            % Bioavailability
P = 1;              % Total drug concentration (normalized to 1)

% Calculate the initial unbound and bound drug concentrations based on the total drug amount
initial_unbound = 1;   % Initial unbound drug concentration
initial_bound = 0;     % Initial bound drug concentration

% Time vector (0 to 24 hours)
tspan = [0 24];

% Define the differential equations for the one-compartment model with plasma protein binding
one_compartment_ode = @(t, y) [
    - k_e * y(1) + k_dissociation * y(2) - k_association * y(1),   % Unbound drug concentration
    k_association * y(1) - k_dissociation * y(2) % Bound drug concentration
];

% Solve the differential equations for the one-compartment model using ode45
[t_one_compartment, concentrations_one_compartment] = ode45(one_compartment_ode, tspan, ...
                                                     [initial_unbound initial_bound]);

% Calculate unbound and bound drug concentrations for the one-compartment model
C_unbound_one_compartment = concentrations_one_compartment(:, 1);
C_bound_one_compartment = concentrations_one_compartment(:, 2);

% Define the differential equations for the two-compartment model without plasma protein binding
k10 = 0.5;          % Rate constant for drug distribution from central to peripheral compartment
k12 = 0.3;          % Rate constant for drug distribution from central to peripheral compartment
k21 = 0.2;          % Rate constant for drug distribution from peripheral to central compartment

two_compartment_ode = @(t, y) [
    -k10 * y(1) + k21 * y(2),                          % Central compartment drug concentration
    k12 * y(1) - k21 * y(2) - k_e * y(2),              % Peripheral compartment drug concentration
    k10 * y(1) - (k12 + k21) * y(3)                    % Drug in transit between central and peripheral compartments
];

% Initial conditions for the two-compartment model
initial_central = 1;   % Initial drug concentration in the central compartment
initial_peripheral = 0; % Initial drug concentration in the peripheral compartment

% Solve the differential equations for the two-compartment model using ode45
[t_two_compartment, concentrations_two_compartment] = ode45(two_compartment_ode, tspan, [initial_central initial_peripheral 0]);

% Calculate drug concentrations for the two-compartment model
C_central = concentrations_two_compartment(:, 1);
C_peripheral = concentrations_two_compartment(:, 2);

% Plot the concentration-time curves for both models on the same plot
figure;

% Plot the one-compartment model with plasma protein binding
plot(t_one_compartment, C_unbound_one_compartment, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Unbound Drug (1-Compartment)');
hold on;
plot(t_one_compartment, C_bound_one_compartment, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Bound Drug (1-Compartment)');

% Plot the two-compartment model without plasma protein binding
plot(t_two_compartment, C_central, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Central Compartment (2-Compartment)');
plot(t_two_compartment, C_peripheral, 'm--', 'LineWidth', 1.5, 'DisplayName', 'Peripheral Compartment (2-Compartment)');

xlabel('Time (hours)');
ylabel('Drug Concentration (mg/L)');
title('Drug Concentrations Comparison between 1-Compartment and 2-Compartment Models');
legend('show');
grid on;