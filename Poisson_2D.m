close all; clc; clear all;
format shortE;
addpath ./forked_MESHND;
addpath ./Morton;
addpath ./Hilbert;
addpath ./Plotting;

q = 5;
n = 2^q
N = n*n;
% n = 2^5 = 32; % Problem size in linear dimension

flag_save_files = 0;

h = 1/(n+1);
h2=(n+1)*(n+1);

K1D = spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);   % 1d Poisson matrix (negative Laplacian)
% subplot(2,3,4), spy(K1D)

I1D = speye(size(K1D));                       % 1d identity matrix
K2D = kron(K1D,I1D)+kron(I1D,K1D);            % 2d Poisson matrix (sparse format)

% % % Natural ordering % % %
p_nat = [1:n*n];		% Natural ordering
A_nat = K2D(p_nat,p_nat);

% % % Algebraic minimum degree % % %
p_amd = amd(K2D);		% AMD ordering
A_amd = K2D(p_amd,p_amd);

% % % Reversed Cuthill-McKee % % %
p_rcm = symrcm(K2D);	% C-M ordering
A_rcm = K2D(p_rcm,p_rcm);

% % % Nested Dissection % % %
G = reshape(1:(n*n*1), n, n, 1)'; 	% Grid
A = -meshsparse(G, 5); 	% 2D stencil - 5 pt stencil
p_ndi = nd2(G); % Get ND permutation
A_ndi = K2D(p_ndi,p_ndi);

% % % Morton ordering (Z) % % %
p_mor = morton(n);
A_mor = K2D(p_mor,p_mor);

% % % Hilbert ordering (H) % % % 
p_hil = hilbert(n)+1; % For 1-based indexing
A_hil = K2D(p_hil,p_hil);


% h = subplot(2,1,1);
% spy(A_nat);
% title(h,'Natural ordering')
% % xlabel(h,'');
% ylabel(h,'');
% set(h,'XTick',[1 256 512 768 1024]);
% set(h,'YTick',[1 256 512 768 1024]);

% h = subplot(2,1,2);
% spy(A_hil);
% title(h,'Hilbert ordering (H)')
% % xlabel(h,'');
% ylabel(h,'');
% set(h,'XTick',[1 256 512 768 1024]);
% set(h,'YTick',[1 256 512 768 1024]);

% return;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Plot permuted matrices
%  Figure size
FigHandle = figure('Position', [100, 100, 650, 650]);

h = subplot(2,3,1);
spy(A_nat);
title(h,'Natural ordering')
xlabel(h,'');
ylabel(h,'');
set(h,'XTick',[]);
set(h,'YTick',[]);
% set(h,'XTick',[1 N*(1/4) N*(1/2) N*(3/4) N]);
% set(h,'YTick',[1 N*(1/4) N*(1/2) N*(3/4) N]);

h = subplot(2,3,2);
spy(A_mor);
title(h,'Morton ordering (Z)')
xlabel(h,'');
ylabel(h,'');
set(h,'XTick',[]);
set(h,'YTick',[]);
% set(h,'XTick',[1 N*(1/4) N*(1/2) N*(3/4) N]);
% set(h,'YTick',[1 N*(1/4) N*(1/2) N*(3/4) N]);

h = subplot(2,3,3);
spy(A_hil);
title(h,'Hilbert ordering (H)')
xlabel(h,'');
ylabel(h,'');
set(h,'XTick',[]);
set(h,'YTick',[]);
% set(h,'XTick',[1 N*(1/4) N*(1/2) N*(3/4) N]);
% set(h,'YTick',[1 N*(1/4) N*(1/2) N*(3/4) N]);

h = subplot(2,3,4);
spy(A_ndi);
title(h,'Nested Dissection')
xlabel(h,'');
ylabel(h,'');
set(h,'XTick',[]);
set(h,'YTick',[]);
% set(h,'XTick',[1 N*(1/4) N*(1/2) N*(3/4) N]);
% set(h,'YTick',[1 N*(1/4) N*(1/2) N*(3/4) N]);

h = subplot(2,3,5);
spy(A_rcm);
title(h,'Cuthill-McKee')
xlabel(h,'');
ylabel(h,'');
set(h,'XTick',[]);
set(h,'YTick',[]);
% set(h,'XTick',[1 N*(1/4) N*(1/2) N*(3/4) N]);
% set(h,'YTick',[1 N*(1/4) N*(1/2) N*(3/4) N]);

h = subplot(2,3,6);
spy(A_amd);
title(h,'Algebraic minimum degree')
xlabel(h,'');
ylabel(h,'');
set(h,'XTick',[]);
set(h,'YTick',[]);
% set(h,'XTick',[1 N*(1/4) N*(1/2) N*(3/4) N]);
% set(h,'YTick',[1 N*(1/4) N*(1/2) N*(3/4) N]);


% Main title
titleString = sprintf('2D Poisson. 5pt stencil. N = %d',N);
p=mtit(titleString,'fontsize',14,'xoff',0,'yoff',.025);

if (flag_save_files ==1)
	disp('Saving matrices to file');
	
	% Save sparse matrices to file
	% save('sA_nat.mat', 'A_nat');
	% save('sA_amd.mat', 'A_amd');
	% save('sA_rcm.mat', 'A_rcm');
	% save('sA_ndi.mat', 'A_ndi');

	% Save Dense matrices into file
	% dA_nat = full(K2D(p_nat,p_nat));
	% dA_amd = full(K2D(p_amd,p_amd));
	% dA_rcm = full(K2D(p_rcm,p_rcm));
	% dA_ndi = full(K2D(p_ndi,p_ndi));
	% save('dA_nat.mat', 'dA_nat', '-ascii');
	% save('dA_amd.mat', 'dA_amd', '-ascii');
	% save('dA_rcm.mat', 'dA_rcm', '-ascii');
	% save('dA_ndi.mat', 'dA_ndi', '-ascii');

	% Save Dense matrices into file

	iA_nat = inv( full( K2D(p_nat,p_nat) ) );
	iA_mor = iA_nat(p_mor,p_mor);
	iA_hil = iA_nat(p_hil,p_hil);
	iA_ndi = iA_nat(p_ndi,p_ndi);
	iA_rcm = iA_nat(p_rcm,p_rcm);
	iA_amd = iA_nat(p_amd,p_amd);

	save('iA_nat.mat', 'iA_nat','-v6');
	save('iA_mor.mat', 'iA_mor','-v6');
	save('iA_hil.mat', 'iA_hil','-v6');
	save('iA_ndi.mat', 'iA_ndi','-v6');
	save('iA_rcm.mat', 'iA_rcm','-v6');
	save('iA_amd.mat', 'iA_amd','-v6');
	disp('End of the program...');
end