
    % Choose which continuous mode n to compare
    mode_index = 3; 
    
    % Compute the continuous mode shape X_n(x)
    x_cont = linspace(0, string_length, 2000);
    X_cont = sin( mode_index*pi*x_cont/string_length );

    % Mass counts to test
    num_masses_list = [10, 30, 100, 300];

    figure('Name','Discrete vs Continuous Mode Shapes','NumberTitle','off');
    hold on;

    colors = lines(length(num_masses_list));

    for k = 1:length(n_list)

        n = num_masses_list(k);
        dx_k = string_length/(n+1);

        % Construct discrete system matrices
        string_params_local = string_params;
        string_params_local.n = n;
        string_params_local.dx = dx_k;

        [M_mat, K_mat] = construct_2nd_order_matrices(string_params_local);

        % Generalized eigenvectors
        [Ur_mat_local, lambda_mat_local] = eig(K_mat, M_mat);

        % Discrete eigenvalue and eigenvector corresponding to mode_index
        lambda_vals = diag(lambda_mat_local);
        [~, idx_sorted] = sort(lambda_vals); % sort eigenvalues
        mode_col = idx_sorted(mode_index);

        Ur_discrete = Ur_mat_local(:, mode_col);

        % Normalize sign (eigenvectors defined only up to ±1)
        % Make its first entry positive so comparisons are consistent
        if Ur_discrete(1) < 0
            Ur_discrete = -Ur_discrete;
        end

        % Rescale discrete mode to match continuous mode magnitude
        X_cont_interp = sin(mode_index*pi*(linspace(0,string_length,n+2)')/string_length);
        X_cont_interp = X_cont_interp(2:end-1);  % internal points only

        Ur_discrete = Ur_discrete * ( norm(X_cont_interp) / norm(Ur_discrete) );

        % Plot discrete approximation
        x_disc = linspace(0, string_length, n+2);
        x_disc = x_disc(2:end-1);  % interior points only

        plot(x_disc, Ur_discrete, 'o-', 'Color', colors(k,:),'LineWidth', 1.5, 'DisplayName', sprintf("Discrete, n = %d masses", n));
    end

    % Plot the continuous mode shape on top
    plot(x_cont, X_cont, 'k--', 'LineWidth', 2.5,'DisplayName', sprintf("Continuous Mode %d", mode_index));

    xlabel("x");
    ylabel("Mode Shape");
    title(sprintf("Mode Shape Comparison for Mode %d", mode_index));
    legend show;
    grid on;
    hold off;