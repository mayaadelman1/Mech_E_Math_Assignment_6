function [omega_vec,Ur_mat] = string_modal_analysis(string_params)


[M_mat,K_mat] = construct_2nd_order_matrices(string_params);

[Ur_mat,lambda_mat] = eig(K_mat,M_mat);


lambda_vec = diag(lambda_mat);
omega_vec  = sqrt(lambda_vec);  


% [omega_vec,idx] = sort(omega_vec);
% Ur_mat = Ur_mat(:,idx);
end
