% Parameters for the two-compartment model with plasma protein binding
k10 = 0.5;          % Rate constant for drug distribution from central to peripheral compartment
k12 = 0.3;          % Rate constant for drug distribution from central to peripheral compartment
k21 = 0.2;          % Rate constant for drug distribution from peripheral to central compartment
k_e = 0.1;          % Rate constant for drug elimination

% Time vector (0 to 24 hours)
tspan = [0 24];

% Initial conditions for drug concentrations and drug-protein complex
initial_conditions = [1 0];  % [Central compartment, Peripheral compartment]

% Time step
dt = 0.1; 

% Number of steps
num_steps = (tspan(2) - tspan(1)) / dt;

% Initialize arrays to store drug concentrations over time
C_central_partition = zeros(num_steps, 1);
C_peripheral_partition = zeros(num_steps, 1);

% Initial conditions
C_central_partition(1) = initial_conditions(1);
C_peripheral_partition(1) = initial_conditions(2);

% Simulate drug concentrations using the partitioning approach
for i = 2:num_steps
    % Update drug concentrations using partitioning
    C_central_partition(i) = C_central_partition(i-1) - k10 * C_central_partition(i-1) * dt + k21 * C_peripheral_partition(i-1) * dt - k12 * C_central_partition(i-1) * dt;
    C_peripheral_partition(i) = C_peripheral_partition(i-1) + k12 * C_central_partition(i-1) * dt - k21 * C_peripheral_partition(i-1) * dt - k_e * C_peripheral_partition(i-1) * dt;
end

% Define the differential equations for the two-compartment model with plasma protein binding
two_compartment_ode_with_binding = @(t, y) [
    -k10 * y(1) + k21 * y(2) - k12 * y(1),                  % Central compartment drug concentration
    k12 * y(1) - k21 * y(2) - k_e * y(2),                    % Peripheral compartment drug concentration
    k10 * y(1) - (k12 + k21) * y(3),           % Drug in transit between central and peripheral compartments
];

% Initial conditions for the two-compartment model
initial_conditions_ode = [1 0 0];  % [Central compartment, Peripheral compartment, Drug in transit]

% Solve the differential equations using ode45
[t_ode45, concentrations_ode45] = ode45(two_compartment_ode_with_binding, tspan, initial_conditions_ode);

% Plot the concentration-time curves for both methods on the same plot
figure;

% Plot the partitioning method
plot(tspan(1) + (0:num_steps-1) * dt, C_central_partition, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Central Compartment (Partitioning)');
hold on;
plot(tspan(1) + (0:num_steps-1) * dt, C_peripheral_partition, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Peripheral Compartment (Partitioning)');

% Plot the ODE45 method
plot(t_ode45, concentrations_ode45(:, 1), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Central Compartment (ODE45)');
plot(t_ode45, concentrations_ode45(:, 2), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Peripheral Compartment (ODE45)');

xlabel('Time (hours)');
ylabel('Drug Concentration');
title('Comparison of Concentration-Time Curves between Partitioning and ODE45 Methods');
legend('show');
grid on;
