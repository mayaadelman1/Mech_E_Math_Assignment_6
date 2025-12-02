
%INPUTS
%t: current time
%V: system state. V = [U;dUdt] where
% U and dUdt are n x 1 column vectors
%string_params: a struct containing the system parameters describing the string
% string_params.n: number of masses
% string_params.M: total mass attached to the string
% string_params.Uf_func: function describing motion of end point
% string_params.dUfdt_func: time derivative of Uf
% string_params.Tf: %tension in string
% string_params.L: %length of string
% string_params.c: %damping coefficient
% string_params.dx: %horizontal spacing between masses
function dVdt = string_rate_func01(t,V,string_params)
    
    n = string_params.n; %number of masses
    M = string_params.M; %total mass attached to the string
    Uf_func = string_params.Uf_func; %function describing motion of end point
    dUfdt_func = string_params.dUfdt_func; %time derivative of Uf
    Tf = string_params.Tf; %tension in string
    L = string_params.L; %length of string
    c = string_params.c; %damping coefficient
    dx = string_params.dx; %horizontal spacing between masses
    
    %unpack state variable
    U = V(1:n);
    dUdt = (V((n+1):(2*n)));

    Uf = Uf_func(t);
    dUfdt = dUfdt_func(t);
    
    % %compute acceleration
    % H = (M/n)*eye(n);
    % I = eye(n);
    % Q = -2*I + circshift(I, [0, 1]) + circshift(I, [0, -1]);
    % Q(1, end) = 0;
    % Q(end, 1) = 0;
    % K = -Tf/dx * Q;
    % B = Tf/dx;


    U_left = [U(2:end);Uf];
    U_right = [0;U(1:end-1)];

    dU_left = [dUdt(2:end);dUfdt];
    dU_right = [0;dUdt(1:end-1)];

    term1 = (Tf/dx)*(U_left-2*U+U_right);
    term2 = c*(dU_left-2*dUdt+dU_right);

    d2Udt2 = (term1+term2)/(M/n);

    %d2Udt2 = (-K * U + B*Uf) \ H;

    % damping = c/dx*dUdt;
    % d2Udt2 = (-K * U + B*Uf + damping) \H;
    % 
    % damping = (c/dx*dUdt)
    % d2Udt2 = ((-K * U + B*Uf + damping) \H)'



    % d2Udt2 = -inv(M)*U
    % K = zeros(n);
    % 
    % for i = 1:n
    %     k(n, n) = -2;
    %     if i > 1 || i < n
    % end
    % %assemble state derivative
    dVdt = [dUdt; d2Udt2];
end