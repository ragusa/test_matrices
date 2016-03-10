close all; clc; clear all;
format shortE;
addpath ./forked_MESHND;
addpath ./Morton;
addpath ./Hilbert;
addpath ./Plotting;

q = 2;
n = 2^q;
N = n

flag_save_files = 1;

K1D = spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);   % 1d Poisson matrix (negative Laplacian)
full(K1D)

% % % Natural ordering % % %
p_nat = [1:n];		% Natural ordering
A_nat = K1D(p_nat,p_nat);

if ( flag_save_files == 1 )
	disp('Saving matrices to file');
	% Save sparse matrices to file
	% save('sA_nat.mat', 'A_nat');

	% Save Dense matrices into file
	iA_nat = inv( full( K1D(p_nat,p_nat) ) );
	save('iA_nat.mat', 'iA_nat','-v6');
	
end

disp('End of the program...');