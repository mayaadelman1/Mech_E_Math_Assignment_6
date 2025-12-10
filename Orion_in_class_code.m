function Orion_in_class_code
    %experiment01();
    %experiment02();
    %experiment03();
    experiment04;
end

function experiment04()   
    string_params = struct();

    n = 5;
    mode_index = 3;
    
    string_params.n = n; %number of masses
    string_params.M = 10; %total mass attached to the string
    string_params.Tf = 2; %tension in string
    string_params.L = 7; %length of string
    string_params.c = 0.0001; %damping coefficient
    string_params.dx = string_params.L/(n+1); %horizontal spacing between masses
    
    rho = string_params.M/string_params.L;
    c = sqrt(string_params.Tf/rho);

    omega_n_spatial = pi*mode_index/string_params.L;
    omega_n = c*omega_n_spatial;
    
    x_list = linspace(0, string_params.L, n+2);
    x_list = x_list(2:end-1);

    x_list_continuous = linspace(0, string_params.L, 1000);
    mode_shape_WE = sin(omega_n_spatial*x_list);
    
    hold on
    plot(x_list_continuous, mode_shape_WE, 'k');
    n = 200;
    string_params.n = n; % number of masses
    string_params.dx = string_params.L/(n+1); %horizontal spacing between masses
    
    n_list = [8, 10, 30, 50];
    
    omega_list = []; 

    for k = 1:length(n_list)
        n = n_list(k);

        [M_mat, K_mat]= construct_2nd_order_matrices(string_params);

        %Use MATLAB to solve the generalized eigenvalue problem
        [Ur_mat,lambda_mat] = eig(K_mat,M_mat);

        mode_shape_LA = Ur_mat([0; Ur_mat(:, n+1-mode_index); 0]);
        mode_shape_LA = 1/max(abs(mode_shape_LA));

        omega_n = sqrt(-lambda_mat(n+1-mode_index, n+1-mode_index));
        omega_list(end+1)
        x_list = linspace(0, string_params.L, n+2)';

        subplot(x_list, mode_shape_LA, 'o-', 'color', 'k', 'linewidth', 2, 'markerfacecolor', 'r', 'markeredgecolor', 'r', 'markersize', 4);
    end
        figure(2);
     
        semilogx(n_list, omega_list, 'o-', 'markerfacecolor', 'r', 'markersize', 3);
        semilogx(n_list, ones(length(n_list), 1)*omega_n_WE, 'b--', 'linewidth', 2);
end


function experiment03()    
    string_params = struct();

    n = 200;
    mode_index = 10;

    string_params.n = n; %number of masses
    string_params.M = 10; %total mass attached to the string
    string_params.Tf = 2; %tension in string
    string_params.L = 7; %length of string
    string_params.c = 0.0001; %damping coefficient
    string_params.dx = string_params.L/(n+1); %horizontal spacing between masses
    
    rho = string_params.M/string_params.L;

    c = sqrt(string_params.Tf/rho);
    
    % [M_mat,K_mat] = construct_2nd_order_matrices(string_params);

    % %Use MATLAB to solve the generalized eigenvalue problem
    % [Ur_mat,lambda_mat] = eig(K_mat,M_mat);
    % 
    % mode_shape_LA = Ur_mat(:, mode_index);

    %omega_n = sqrt(-lmabda_mat(mode_index, mode_index));

    omega_n_spatial = pi*mode_index/string_params.L;
    omega_n = c*omega_n_spatial;
    
    x_list = linspace(0, string_params.L, n+2);
    x_list = x_list(2:end-1);

    mode_shape_WE = sin(omega_n_spatial*x_list);
    

    omega = omega_n;
    A = 3;

    string_params.Uf_func = @(t_in) A*sin(omega*t_in); %function describing motion of end point
    string_params.dUfdt_func = @(t_in) -omega*A*cos(omega*t_in); %time derivative of Uf
    string_params.Uf_amplitude = A;

    U0 = zeros(n, 1);
    dUdt0 = zeros(n,1);

    V0 = [U0, dUdt0];

    tlist = linspace(0, 20*(2*pi)/omega, 1000+1);
    my_rate_func = @(t_in, V_in) string_rate_func01(t_in, V_in, string_params);
    
    [tlist, Vresult] = ode45(my_rate_func, tlist, V0);
    %animate_string_with_mode_shape(tlist, Vresult,string_params,mode_shape_WE);
    animate_string_travelling_wave(tlist, Vresult, string_params, x_centroid0, c);
end

function animate_string_travelling_wave(tlist, Vresult, string_params, x_centroid0, c)
    n = stringparams.n;
    L = stringparams.L;

    x_list = linspace(0, L, n+2);

    % string_plot = plot(0, 0, 'o-', 'color', 'k', 'linewidth', 2, 'markerfacecolor', 'r', 'markeredgecolor', 'r', 'markersize', 4);
    % 
    % travelling_wave_plot = plot(0, 0, 'o-', 'color', 'k', 'linewidth', 2, 'markerfacecolor', 'r', 'markeredgecolor', 'r', 'markersize', 4);
    % 
    % travelling_wave_padded = [0;mode_shape_LA;0];
    % 
    % scale_factor = max*abs(maxU);
    % 
    % max(abs(mode_shape_padded));
    % set(mode_shape_plot, 'xdata', x-list,' ydata', mode_shape_padded);

    
    maxU = max(max(Vresult(:, 1:n)));
    minU = min(min(Vresult(:, 1:n)));

    h = max([abs(maxU), abs(minU), string_params.Uf_amplitude]);
    
    axis([0, L, -1.1*h, 1.1*h]);
    hold on

    string_plot = plot(0, 0, 'o-', 'color', 'k', 'linewidth', 2, 'markerfacecolor', 'r', 'markeredgecolor', 'r', 'markersize', 4);
    centroid_plot = plot(0,0, 'b', 'linewidth', 2);

    xlabel('x');
    ylabel('U(t,x)');
    title('Plot of vibrating string');
   

    for k = 1:length(tlist)

        % %short blurb showing how to find x-coord of tracking line
        % %x = x-coord of tracking line, t = current time
        % %c = wave speed, w = pulse width (in time), L = string length
        % x_centroid = c*t+x_centroid0;
        % x_centroid = mod(x_centroid, 2*L);
        % if x_centroid > L
        %     x_centroid = 2*L - x_centroid;
        % end
        % 
        % set(centroid_plot, 'xdata', [x_centroid, x_centroid], 'y_data', [-1.1*h, 1.1*h]);
        t = tlist(k);
        U = Vresult(k, 1:n);
        Uf = string_params.Uf_func(t);
        
        U_padded = [0, U, Uf];
        x_list = linspace(0, L, n+2);
        
        set(string_plot, 'xdata', x_list, 'ydata', U_padded);
        drawnow;
    end
end

function experiment02()
    omega = 5;
    A = 3;
    
    string_params = struct();

    n = 200;
    mode_index = 4;

    string_params.n = 5; %number of masses
    string_params.M = 10; %total mass attached to the string
    string_params.Tf = 2; %tension in string
    string_params.L = 7; %length of string
    string_params.c = 0.0001; %damping coefficient
    string_params.dx = string_params.L/(n+1); %horizontal spacing between masses


    string_params.Uf_func = @(t_in) 0; %function describing motion of end point
    string_params.dUfdt_func = @(t_in) 0; %time derivative of Uf
    string_params.Uf_amplitude = 0;

    rho = string_params.M/string_params.L;

    c = sqrt(string_params.Tf/rho);
    
    x_list = linspace(0, string_params.L, n+2)';
    x_list = x_list(2:end-1);
    
    w_pulse = string_params.L/5;
    h_pulse = 7;

    x_shift = 0;
    U0 = bspline_pulse(x_list-x_shift, w_pulse, h_pulse);
    dUdt0 = -c*bspline_pulse_derivative(x_list-x_shift, w_pulse, h_pulse);
    
    figure(1)
    hold on
    plot(x_list, U0, 'o-', 'color', 'r');
    plot(x_list, dUdt, 'o-', 'color', 'b');

    figure(2)
    V0 = [U0, dUdt0];

    tlist = linspace(0, 0.0001*3*string_params.L/c, 5000+1);

    my_rate_func = @(t_in, V_in) string_rate_func01(t_in, V_in, string_params);

    [tlist, Vresult] = ode45(my_rate_func, tlist, V0);

    animate_string_travelling_wave(tlist, Vresult,string_params,x_centroid0, c);
end

% %triangle pulse function
% %INPUTS:
% %t: current time
% %w: width of pulse (starts at t=0, ends at t=w)
% %h: height of pulse
% %OUTPUTS:
% %res: pulse evaluated at t
% function res = triangle_pulse(t,w,h)
% t = t*(2/w);
% res = 1-min(1*abs(t-1),1);
% res = h*res;
% end
% 
% %triangle pulse function (derivative)
% %INPUTS:
% %t: current time
% %w: width of pulse (starts at t=0, ends at t=w)
% %h: height of pulse
% %OUTPUTS:
% %res: derivative of pulse evaluated at t
% function res = triangle_pulse_derivative(t,w,h)
% t = t*(2/w);
% res = -sign(t-1).*(abs(t-1)<1);
% res = (2*h/w)*res;
% end

%b-spline pulse function
%INPUTS:
%t: current time
%w: width of pulse (starts at t=0, ends at t=h)
%h: height of pulse
%OUTPUTS:
%res: pulse evaluated at t
function res = b_spline_pulse(t,w,h)
t = 4*t/w;
b3 = (0<=t).*(t<1).*(t.^3)/4;
t = t-1;
b2 = (0<=t).*(t<1).*(-3*t.^3+3*t.^2+3*t+1)/4;
t = t-1;
b1 = (0<=t).*(t<1).*(3*t.^3-6*t.^2+4)/4;
t = t-1;
b0 = (0<=t).*(t<1).*(-t.^3+3*t.^2-3*t+1)/4;
res = h*(b0+b1+b2+b3);
end
%b-spline pulse function (derivative)
%INPUTS:
%t: current time
%w: width of pulse (starts at t=0, ends at t=h)
%h: height of pulse
%OUTPUTS:
%res: derivative of pulse evaluated at t

function res = b_spline_pulse_derivative(t,w,h)
t = 4*t/w;
b3 = (0<=t).*(t<1).*(3*t.^2)/4;
t = t-1;
b2 = (0<=t).*(t<1).*(-9*t.^2+6*t+3)/4;
t = t-1;
b1 = (0<=t).*(t<1).*(9*t.^2-12*t)/4;
t = t-1;
b0 = (0<=t).*(t<1).*(-3*t.^2+6*t-3)/4;
res = (4*h/w)*(b0+b1+b2+b3);
end

function experiment01()
    omega = 5;
    A = 3;
    
    string_params = struct();

    n = 5;
    mode_index = 4;

    string_params.n = 5; %number of masses
    string_params.M = 10; %total mass attached to the string
    string_params.Tf = 2; %tension in string
    string_params.L = 7; %length of string
    string_params.c = 0.0001; %damping coefficient
    string_params.dx = string_params.L/(n+1); %horizontal spacing between masses

    [M_mat,K_mat] = construct_2nd_order_matrices(string_params);

    %Use MATLAB to solve the generalized eigenvalue problem
    [Ur_mat,lambda_mat] = eig(K_mat,M_mat);

    mode_shape_LA = Ur_mat(:, mode_index);

    omega_n = sqrt(-lmabda_mat(mode_index, mode_index));

    omega = omega_n;
    A = 3;

    string_params.Uf_func = @(t_in) A*sin(omega*t_in); %function describing motion of end point
    string_params.dUfdt_func = @(t_in) -omega*A*cos(omega*t_in); %time derivative of Uf
    string_params.Uf_amplitude = A;

    U0 = zeros(n, 1);
    dUdt0 = zeros(n,1);

    V0 = [U0, dUdt0];

    tlist = linspace(0, 20*(2*pi)/omega, 1000+1);
    my_rate_func = @(t_in, V_in) string_rate_func01(t_in, V_in, string_params);
    
    [tlist, Vresult] = ode45(my_rate_func, tlist, V0);
    animate_string_with_mode_shape(tlist, Vresult,string_params,mode_shape_LA);
end


function animate_string_with_mode_shape(tlist, Vresult, string_params)
    n = stringparams.n;
    L = stringparams.L;
    x_list = linspace(0, L, n+2);

    string_plot = plot(0, 0, 'o-', 'color', 'k', 'linewidth', 2, 'markerfacecolor', 'r', 'markeredgecolor', 'r', 'markersize', 4);
    
    mode_shape_plot = plot(0, 0, 'o-', 'color', 'k', 'linewidth', 2, 'markerfacecolor', 'r', 'markeredgecolor', 'r', 'markersize', 4);
    
    mode_shape_padded = [0;mode_shape_LA;0];

    scale_factor = max*abs(maxU);

    max(abs(mode_shape_padded));
    set(mode_shape_plot, 'xdata', x-list,' ydata', mode_shape_padded);

    maxU = max(max(Vresult(:, 1:n)));
    minU = min(min(Vresult(:, 1:n)));

    h = max(abs(maxU), abs(minU));
    
    axis([0, L, -1.1*h, 1.1*h]);
    hold on
    xlabel('x');
    ylabel('U(t,x)');
    title('Plot of vibrating string');

    for k = 1:length(tlist)
        t = tlist(k);
        U = Vresult(k, 1:n);
        Uf = string_params.Uf_func(t);
        
        U_padded = [0, U, Uf];
        x_list = linspace(0, L, n+2);
        
        set(string_plot, 'xdata', x_list, 'ydata', U_padded);
        drawnow;
    end
end


%build the mass and stiffness matrices that describe the 2nd order system.
%INPUTS
%string_params: a struct containing the system parameters describing the string
% string_params.n: number of masses
% string_params.M: total mass attached to the string
% string_params.Uf_func: function describing motion of end point
% string_params.dUfdt_func: time derivative of Uf
% string_params.Tf: %tension in string
% string_params.L: %length of string
% string_params.c: %damping coefficient
% string_params.dx: %horizontal spacing between masses
%OUTPUTS
%M_mat: the n x n mass (inertia) matrix
%K_mat: the n x n stiffness matrix
function [M_mat,K_mat] = construct_2nd_order_matrices(string_params)
    Tf = string_params.Tf; %tension in string
    L = string_params.L; %length of string
    c = string_params.c; %damping coefficient
    dx = string_params.dx; %horizontal spacing between masses
    
    % construct the nxn discrete Laplacian matrix
    n = string_params.n;
    I_n = eye(n);
    M_left = circshift(I_n, [0, -1]);
    M_right = circshift(I_n, [0, 1]);
    my_Laplacian = M_left-2*I_n+M_right;

    my_Laplacian(1,n) = 0;
    my_Laplacian(n,1) = 0;
    
    K_mat = (Tf/dx)*my_Laplacian;
    M_mat = (M/n)*I_n;

end

