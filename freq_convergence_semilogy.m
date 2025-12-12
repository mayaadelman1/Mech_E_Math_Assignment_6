function [N_list,omega_disc,omega_cont,err] = freq_convergence_semilogy(string_params,k_mode,N_list,use_relative)
if nargin < 2 || isempty(k_mode), k_mode = 1; end
if nargin < 3 || isempty(N_list), N_list = [5 10 20 40 80 160 320]; end
if nargin < 4 || isempty(use_relative), use_relative = false; end

rho = string_params.M/string_params.L;
c = sqrt(string_params.Tf/rho);
omega_cont = c*pi*k_mode/string_params.L;

omega_disc = zeros(size(N_list));
for i = 1:numel(N_list)
    N = N_list(i);
    params = string_params;
    params.n = N;
    params.dx = params.L/(N+1);
    [M_mat,K_mat] = construct_2nd_order_matrices(params);
    [~,lambda_mat] = eig(K_mat,M_mat);
    omega = sqrt(abs(diag(lambda_mat)));
    omega = sort(omega);
    omega_disc(i) = omega(k_mode);
end

err = abs(omega_disc - omega_cont);
if use_relative
    err = err/omega_cont;
end

figure;
semilogy(N_list,err,'o-','LineWidth',1.5);
grid on;
xlabel('Number of masses, N');
if use_relative
    ylabel('Relative error |ω_{disc} − ω_{cont}|/ω_{cont}');
else
    ylabel('Absolute error |ω_{disc} − ω_{cont}|');
end
title(sprintf('Convergence of Mode %d Resonant Frequency',k_mode));
end
