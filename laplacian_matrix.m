function my_Laplacian = laplacian_matrix(n)
    %construct the nxn discrete laplacian matrix
    %H = (M/n)*eye(n);
    I_n = eye(n); % build the nxn identity matrix
    my_Laplacian = -2*I_n + circshift(I_n, [0, 1]) + circshift(I_n, [0, -1]);
    my_Laplacian(1,end) = my_Laplacian(1,end)-1; %delete unwanted 1 in top right corner
    my_Laplacian(end,1) = my_Laplacian(end,1)-1; %delete unwanted 1 in bottom right corner
    %K = -Tf/dx * Q;
    %B = Tf/dx;
end
