close all; clc; clear all;
addpath ./forked_MESHND;

n = 32; % Problem size in linear dimension
flag_save_files = 0;

h = 1/(n+1);
h2=(n+1)*(n+1);

K1D = spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);   % 1d Poisson matrix (negative Laplacian)
% subplot(2,3,4), spy(K1D)

I1D = speye(size(K1D));                       % 1d identity matrix
K2D = kron(K1D,I1D)+kron(I1D,K1D);            % 2d Poisson matrix (sparse format)

% Computes permutation ordering
p_nat = [1:1024];		% Natural ordering
p_amd = amd(K2D);		% AMD ordering
p_rcm = symrcm(K2D);	% C-M ordering

% NDi ordering
G = reshape(1:(n*n*1), n, n, 1)'; 	% Grid
A = -meshsparse(G, 5); 	% 2D stencil - 5 pt stencil
p_ndi = nd2(G); % Get ND permutation

% Permute matrices
A_nat = K2D(p_nat,p_nat);
A_amd = K2D(p_amd,p_amd);
A_rcm = K2D(p_rcm,p_rcm);
A_ndi = K2D(p_ndi,p_ndi);

% Plot permuted matrices
%  Figure size
FigHandle = figure('Position', [100, 100, 650, 650]);

h = subplot(2,2,1);
spy(A_nat);
title(h,'Natural ordering')
xlabel(h,'');
ylabel(h,'');
set(h,'XTick',[1 256 512 768 1024]);
set(h,'YTick',[1 256 512 768 1024]);

h = subplot(2,2,2);
spy(A_amd);
title(h,'Algebraic minimum degree (AMD)')
xlabel(h,'');
ylabel(h,'');
set(h,'XTick',[1 256 512 768 1024]);
set(h,'YTick',[1 256 512 768 1024]);

h = subplot(2,2,3);
spy(A_rcm);
title(h,'Reversed Cuthill-McKee (rCM)')
xlabel(h,'');
ylabel(h,'');
set(h,'XTick',[1 256 512 768 1024]);
set(h,'YTick',[1 256 512 768 1024]);

h = subplot(2,2,4);
spy(A_ndi);
title(h,'Nested Dissection (ND)')
xlabel(h,'');
ylabel(h,'');
set(h,'XTick',[1 256 512 768 1024]);
set(h,'YTick',[1 256 512 768 1024]);

% Main title
p=mtit('2D Poisson. 5pt stencil. N = 1,024',...
     'fontsize',14,...
     'xoff',0,'yoff',.025);

if (flag_save_files ==1)
	disp('Saving matrices to file');
	
	% Save sparse matrices to file
	save('sA_nat.mat', 'A_nat');
	save('sA_amd.mat', 'A_amd');
	save('sA_rcm.mat', 'A_rcm');
	save('sA_ndi.mat', 'A_ndi');

	% Save Dense matrices into file
	dA_nat = full(K2D(p_nat,p_nat));
	dA_amd = full(K2D(p_amd,p_amd));
	dA_rcm = full(K2D(p_rcm,p_rcm));
	dA_ndi = full(K2D(p_ndi,p_ndi));
	save('dA_nat.mat', 'dA_nat', '-ascii');
	save('dA_amd.mat', 'dA_amd', '-ascii');
	save('dA_rcm.mat', 'dA_rcm', '-ascii');
	save('dA_ndi.mat', 'dA_ndi', '-ascii');

	% Save Dense matrices into file
	iA_nat = inv( full( K2D(p_nat,p_nat) ) );
	iA_amd = inv( full( K2D(p_amd,p_amd) ) );
	iA_rcm = inv( full( K2D(p_rcm,p_rcm) ) );
	iA_ndi = inv( full( K2D(p_ndi,p_ndi) ) );
	save('iA_nat.mat', 'iA_nat', '-ascii');
	save('iA_amd.mat', 'iA_amd', '-ascii');
	save('iA_rcm.mat', 'iA_rcm', '-ascii');
	save('iA_ndi.mat', 'iA_ndi' ,'-ascii');
end