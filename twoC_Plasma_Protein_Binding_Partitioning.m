function drug_model_with_binding()

    % Parameters and initial conditions
    k10 = 0.1;   % Elimination rate constant from central compartment
    k12 = 0.2;   % Rate constant for drug distribution from central to peripheral compartment
    k21 = 0.15;  % Rate constant for drug distribution from peripheral to central compartment
    kon = 0.01;  % Association rate constant
    koff = 0.05; % Dissociation rate constant
    
    Vc = 1;      % Volume of central compartment
    Vp = 2;      % Volume of peripheral compartment
    P = 10;      % Concentration of plasma protein
    
    % Initial drug concentrations
    Dc_initial = 100;
    Dp_initial = 0;
    DP_initial = 0;
    
    % Time span for simulation
    tspan = [0 10]; % Simulation time span
    
    % Solve the differential equations using ode45
    [t, y] = ode45(@ode_fun, tspan, [Dc_initial; Dp_initial; DP_initial]);
    
    % Extract drug concentrations from the solution
    Dc = y(:, 1);
    Dp = y(:, 2);
    DP = y(:, 3);
    
    % Plot drug concentrations
    figure;
    plot(t, Dc, 'r', t, Dp, 'b', t, DP, 'g');
    xlabel('Time');
    ylabel('Drug Concentration');
    legend('Central Compartment (Free)', 'Peripheral Compartment', 'Bound Drug');
    title('Drug Concentrations in the Two-Compartment Model with Binding');
    
    % Define the system of ODEs
    function dydt = ode_fun(t, y)
        Dc = y(1); % Drug concentration in the central compartment
        Dp = y(2); % Drug concentration in the peripheral compartment
        DP = y(3); % Drug concentration in the bound drug compartment
        
        % Plasma protein binding equation
        f_B = DP / (DP + Dc);
        
        % Rate equations for the three compartments
        dydt = [
            -k10/Vc * Dc - k12/Vc * Dc + k21/Vp * Dp - kon * Dc * P + koff * DP; % Rate of change in central compartment
            k12/Vc * Dc - k21/Vp * Dp; % Rate of change in peripheral compartment
            kon * Dc * P - koff * DP; % Rate of change in bound drug compartment
        ];
    end

end
