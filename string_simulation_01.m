function string_simulation_01()
    clc;
    num_masses = 250;
    total_mass = num_masses;
    tension_force = 2;
    string_length = 0.5;
    damping_coeff = 0.001;
    dx = string_length/(num_masses+1);
    amplitude_Uf = 0.5;
    omega_Uf = 0.9922;

    %list of x points (including the two endpoints)
    xlist = linspace(0,string_length,num_masses+2);
    Uf_func = @(t_in) amplitude_Uf*cos(omega_Uf*t_in);
    dUfdt_func = @(t_in) -omega_Uf*amplitude_Uf*sin(omega_Uf*t_in);


    w = string_length/5;
    h = 2;
   
    
    %generate the struct
    string_params = struct();
    string_params.n = num_masses;
    string_params.M = total_mass;
    string_params.Tf = tension_force;
    string_params.L = string_length;
    string_params.c = damping_coeff;
    string_params.dx = dx;
    
    c_wave = sqrt(tension_force * string_length / total_mass); % wave at speed c for wave equation
    % 
    % Uf_func = @(t_in) triangle_pulse(t_in, w, h);
    % dUfdt_func = @(t_in) c_wave*triangle_pulse_derivative(t_in, w, h);
    
    Uf_func = @(t_in) b_spline_pulse(t_in, w/c_wave, h);
    dUfdt_func = @(t_in) c_wave*b_spline_pulse_derivative(t_in, w/c_wave, h);

    % Uf_func = @(t_in) 0*t_in;
    % dUfdt_func = @(t_in) 0*t_in;

    string_params.Uf_func = Uf_func;
    string_params.dUfdt_func = dUfdt_func;
    
    % run modal analysis
    [omega_vec, Ur_mat] = string_modal_analysis(string_params);
    omega_vec(1)


    %load string_params into rate function
    my_rate_func = @(t_in,V_in) string_rate_func01(t_in,V_in,string_params);



    %initial conditions
    U0 = zeros(num_masses,1);
    dUdt0 = zeros(num_masses,1);

    % x_temp = linspace(0,string_params.L,num_masses+2)';
    % x_temp = x_temp(2:end-1);

    % U0 = b_spline_pulse(x_temp, w, h);
    % dUdt0 =  -c_wave*b_spline_pulse_derivative(x_temp, w, h);
    V0 = [U0;dUdt0];
    size(V0)
    % tspan = [0 25];

    tspan = linspace(0,3*string_params.L/c_wave,5001);
    %run the integration
    h_ref = 0.002; BT_struct = get_BT("Dormand Prince");

    [tlist,Vlist] = ode45(my_rate_func,tspan,V0);
    % [tlist,Vlist] = explicit_RK_fixed_step_integration(my_rate_func,tspan,V0, h_ref, BT_struct);
    %disp(xlist)
    %disp(Vlist)
    

    % Mode Shape Plotting

    mode_points = Ur_mat(:,1);
    mode_points = [0; mode_points; 0];

    fig1 =figure('Name','Mode Shape','NumberTitle','off'); hold on
    ylim([-1 1])

    masses = scatter(xlist, mode_points, "filled");
    string = plot(xlist, mode_points, 'LineWidth', 2);
    hold off;

    % Animation
    %fig2 = figure('Name','Animation','NumberTitle','off');
    hold on
    ylim([-5 5])
    tracking_line = xline(0, 'r', 'LineWidth', 2);

    points = Vlist(:, 1:num_masses);
    points = [zeros(length(points), 1) points Uf_func(tlist)];

    %masses = scatter(xlist, points(1,:), "filled");
    string = plot(xlist, points(1,:), 'LineWidth', 2);

    tdiff = [0; diff(tlist)];
    
    for i = 1:length(tlist)
        %set(masses, 'XData', xlist, 'YData', points(i,:));
        %set(string, 'XData', xlist, 'YData', points(i,:));
        
        % tracking line position
        x_pos = string_length + c_wave * tlist(i) - w/2;       % pulse moving to the right
        x_pos = mod(x_pos, 2*string_length);

        if x_pos > string_length
            x_pos = 2*string_length - x_pos;    
        end
        
        set(tracking_line, 'Value', x_pos);

        %drawnow;
        title("Traveling Wave Animation");
        %disp(tlist(i));
        %pause(tdiff(i));
    end


    % MODAL ANALYSIS

    % Mode shapes and frequencies
    rho = total_mass / string_length;          % linear density
    num_modes = 6;                               % number of continuous modes to plot

    x_cont = linspace(0, string_length, 2000); 

    figure('Name','Continuous Mode Shapes','NumberTitle','off'); 
    tiledlayout(3,2);

    for n = 1:num_modes
    
        % continuous resonant frequency
        omega_n = c_wave * (pi*n) / string_length; 

        % continuous mode shape
        Xn = sin( (pi*n/string_length) * x_cont ); 

        nexttile
        plot(x_cont, Xn, 'LineWidth', 2);
        title(sprintf("Continuous Mode %d  (ω_n = %.3f rad/s)", n, omega_n));
        xlabel("x"); ylabel("X_n(x)");
    end

    sgtitle('Continuous Wave Equation Mode Shapes and Frequencies');

 
   % Mode shape vs. discrete approximations comparisons
   
   % PSEUDOCODE FOR THIS SECTION:
   % 1) make num masses list (5-7 ish values)
   % 2) Set n and .dx string params to what you want them to be (experiment04 in Orion's code)
   % 3) Construct M and K mat from string_modal_analysis
   % 4) Compute Ur mat and lamda mat from that
   % 5) Take your mode shape code (this line:  Xn = sin( (pi*n/string_length) * x_cont ));. That will be one thing you're plotting.
   % 6) Add a zero to left and right of Ur_mat to include the endpoints
   % 7) Rescale mode shape la
        %Ex: V = V/max(abs(V)); for any vector v
   
   % 8) Compute corresponding frequency (omega) using lambda mat from eig function in string_modal_analysis. Store result in list. 
   % 9) Generate x-coordinates from masses based on your num masses list. 

      % 8) and 9) is basically this code from Orion's in class code:
            % omega_n_spatial = pi*mode_index/string_params.L;
            % omega_n = c*omega_n_spatial;
            % 
            % 
            % 
            % x_list_continuous = linspace(0, string_params.L, 1000);
            % mode_shape_WE = sin(omega_n_spatial*x_list);

    n = 5;
    mode_index = 3;

    string_params.n = n; % number of masses
    string_params.dx = string_params.L/(n+1); %horizontal spacing between masses
    
    num_masses = [8, 10, 30, 50];

    for k = 1:length(num_masses)
        n = num_masses(k);

        [M_mat, K_mat]= construct_2nd_order_matrices(string_params);

        %Use MATLAB to solve the generalized eigenvalue problem
        [Ur_mat,lambda_mat] = eig(K_mat,M_mat);

        mode_shape_LA = Ur_mat([0; Ur_mat(:, n+1-mode_index); 0]);
        mode_shape_LA = 1/max(abs(mode_shape_LA));

        omega_n_spatial = pi*mode_index/string_params.L;
        omega_n = c*omega_n_spatial;
    
        x_list = linspace(0, string_params.L, n+2);
        x_list = x_list(2:end-1);

        x_list_continuous = linspace(0, string_params.L, 1000);
        mode_shape_WE = sin(omega_n_spatial*x_list);

        subplot(x_list, mode_shape_LA, 'o-', 'color', 'k', 'linewidth', 2, 'markerfacecolor', 'r', 'markeredgecolor', 'r', 'markersize', 4);
    end
end

    