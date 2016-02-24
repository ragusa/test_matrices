
% Change n, to increase matrix size
n = 4;

h = 1/(n+1);
h2=(n+1)*(n+1);

K1D = spdiags(ones(n,1)*[-1 2 -1],-1:1,n,n);  % 1d Poisson matrix
% subplot(2,3,4), spy(K1D)

I1D = speye(size(K1D));                       % 1d identity matrix
K2D = kron(K1D,I1D)+kron(I1D,K1D);            % 2d Poisson matrix (sparse format)
spy(K2D);

f2D = h^2*ones(n^2,1);                        % 2d right hand side
u2D = K2D\f2D;                                % exact solution