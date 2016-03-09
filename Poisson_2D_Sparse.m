close all; clc; clear all;
format shortE;
addpath ./forked_MESHND;
addpath ./Morton;
addpath ./Hilbert;
addpath ./Plotting;

q = 9;
n = 2^q
N = n*n
% n = 2^9 = 512; % Problem size in linear dimension

flag_save_files = 1;

h = 1/(n+1);
h2=(n+1)*(n+1);

K1D = spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);   % 1d Poisson matrix (negative Laplacian)
% subplot(2,3,4), spy(K1D)

I1D = speye(size(K1D));                       % 1d identity matrix
K2D = kron(K1D,I1D)+kron(I1D,K1D);            % 2d Poisson matrix (sparse format)
% spy(K2D)
% % % % Natural ordering % % %
% p_nat = [1:n*n];		% Natural ordering
% A_nat = K2D(p_nat,p_nat);

% % % % Algebraic minimum degree % % %
% p_amd = amd(K2D);		% AMD ordering
% A_amd = K2D(p_amd,p_amd);

% % % % Reversed Cuthill-McKee % % %
% p_rcm = symrcm(K2D);	% C-M ordering
% A_rcm = K2D(p_rcm,p_rcm);

% % % % Nested Dissection % % %
% G = reshape(1:(n*n*1), n, n, 1)'; 	% Grid
% A = -meshsparse(G, 5); 	% 2D stencil - 5 pt stencil
% p_ndi = nd2(G); % Get ND permutation
% A_ndi = K2D(p_ndi,p_ndi);

% % % % Morton ordering (Z) % % %
% p_mor = morton(n);
% A_mor = K2D(p_mor,p_mor);

% % % % Hilbert ordering (H) % % % 
% p_hil = hilbert(n)+1; % For 1-based indexing
% A_hil = K2D(p_hil,p_hil);

save('K2D_262k.mat', 'K2D','-v6');
