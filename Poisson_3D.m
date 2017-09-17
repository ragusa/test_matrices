clc; clear all;

n = 2;
h = 1/(n+1);
h2=(n+1)*(n+1);

K1D = spdiags(ones(n,1)*[-1 2 -1],-1:1,n,n);  % 1d Poisson matrix

I1D = speye(size(K1D));                       % 1d identity matrix
K2D = kron(K1D,I1D)+kron(I1D,K1D);            % 2d Poisson matrix

I2D = speye(size(K2D));                       % 2d identity matrix
K3D = kron(K2D,I1D)+kron(I2D,K1D);            % 3d Poisson matrix

return;

% spy(K3D);

FILENAME_MAT = strcat('K3D_n',int2str(n),'.mat')
% save(FILENAME_MAT, 'K3D');

FILENAME_MTX = strcat('K3D_n',int2str(n),'.mtx')
mmwrite(FILENAME_MTX, K3D, 'GC: Saving K3D matrix')
