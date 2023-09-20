% One-Compartment Model for Acetaminophen (APAP)
CL = 18.7;  % Clearance in l/h 
Vd = 30.9;  % Distribution Volume in l 
k = CL/Vd; % Elimination rate constant 
t = 0:0.1:24;  % Time points from 0 to 24 hours
C0 = 300;     % Initial drug amount 

% Solve the differential equation
C_one_compartment = C0 * exp(-k * t);

% Plot the drug amount over time
figure;
plot(t, C_one_compartment, 'b');
xlabel('Time (hours)');
ylabel('APAP amount (mg)');
title('One-Compartment Model for Acetaminophen (APAP)');
