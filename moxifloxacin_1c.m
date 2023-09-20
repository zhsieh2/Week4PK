% One-Compartment Model for Moxifloxacin
% Define parameters
k = 0.25; % Elimination rate constant 
% Define the time span for the simulation
tspan = [0 24]; % Start at time 0 and end at time 24
% Initial drug amount
C0 = 300; % Initial amount (mg)
% Solve the differential equation using ode45 and the oneCompModel function
[t,C] = ode45(@(t,C) (-k)*C,  tspan, C0);
% Plot the drug amount vs. time
plot(t, C);
xlabel('Time (hours)');
ylabel('Drug Amount (mg)');
title('One-Compartment Pharmacokinetic Model for Moxifloxacin');
grid on;