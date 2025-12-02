function laplacian_matrix
    H = (M/n)*eye(n);
    I_n = eye(n);
    Q = -2*I + circshift(I, [0, 1]) + circshift(I, [0, -1]);
    Q(1, end) = 0;
    Q(end, 1) = 0;
    K = -Tf/dx * Q;
    B = Tf/dx;
    %construct the nxn discrete laplacian matrix
    n = % your code here
    I_n = eye(n); % build the nxn identity matrix
    my_Laplacian = %your code here
    %doing this as subtracting deals with the weird edge cases
    %that occur when n=1 and n=2
    my_Laplacian(1,end) = my_Laplacian(1,end)-1; %delete unwanted 1 in top right corner
    my_Laplacian(end,1) = my_Laplacian(end,1)-1; %delete unwanted 1 in bottom right corner
end
