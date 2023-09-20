function drug_model()

    % Parameters and initial conditions
    k10 = 0.1;   % Elimination rate constant from central compartment
    k12 = 0.2;   % Rate constant for drug distribution from central to peripheral compartment
    k21 = 0.15;  % Rate constant for drug distribution from peripheral to central compartment
    kon = 0.01;  % Association rate constant of plasma protein binding
    koff = 0.05; % Dissociation rate constant of plasma protein binding
    
    Vc = 1;      % Volume of central compartment
    Vp = 2;      % Volume of peripheral compartment
    P = 10;      % Concentration of plasma protein
    
    % Initial drug concentrations
    Dc_initial = 100;
    Dp_initial = 0;
    
    % Time span for simulation
    tspan = [0 10]; % Simulation time span
    
    % Solve the differential equations using ode45
    [t, y] = ode45(@ode_fun, tspan, [Dc_initial; Dp_initial]);
    
    % Extract drug concentrations from the solution
    Dc = y(:, 1);
    Dp = y(:, 2);
    
    % Plot drug concentrations
    figure();
    plot(t, Dc, 'r', t, Dp, 'b');
    xlabel('Time');
    ylabel('Drug Concentration');
    legend('Central Compartment', 'Peripheral Compartment');
    title('Drug Concentrations in the Two-Compartment Model with Effect of Plasma Protein Binding');
    
    % Define the system of ODEs
    function dydt = ode_fun(t, y)
        Dc = y(1); % Drug concentration in the central compartment
        Dp = y(2); % Drug concentration in the peripheral compartment
        
        % Association and dissociation equations
        DP = kon * Dc * P / (koff + kon); % Bound drug concentration
        
        % Rate equations for the two compartments
        dydt = [
            -k10/Vc * Dc - k12/Vc * Dc + k21/Vp * Dp; % Rate of change in central compartment
            k12/Vc * Dc - k21/Vp * Dp;                  % Rate of change in peripheral compartment
        ];
    end

end
