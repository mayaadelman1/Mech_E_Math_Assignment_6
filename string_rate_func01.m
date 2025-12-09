
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
   
    U_left = [U(2:end);Uf];
    U_right = [0;U(1:end-1)];

    dU_left = [dUdt(2:end);dUfdt];
    dU_right = [0;dUdt(1:end-1)];

    term1 = (Tf/dx)*(U_left-2*U+U_right);
    term2 = (c/dx)*(dU_left-2*dUdt+dU_right);

    d2Udt2 = (term1+term2)/(M/n);

    % computing boundary condition special cases
    
    % % first point:
    % term1 = (Tf/dx)*(0 - 2*U(1) + U(2));
    % term2 = (c/dx)*(0 - 2*dUdt(1) + dUdt(2));
    % d2Udt2(1) = (term1+term2)/(M/n);
    % 
    % % last point
    % term1 = (Tf/dx) * (U(n-1) - 2*U(n) + Uf);
    % term2 = (c/dx) * (dUdt(n-1) - 2*dUdt(n) + dUfdt);
    % d2Udt2(n) = (term1+term2)/(M/n);

    %d2Udt2 = (-K * U + B*Uf) \ H;

    % damping = c/dx*dUdt;
    % d2Udt2 = (-K * U + B*Uf + damping) \H;
    % 
    % damping = (c/dx*dUdt)
    % d2Udt2 = ((-K * U + B*Uf + damping) \H)'

    dVdt = [dUdt; d2Udt2];
end