close all; clc; clear all;
format shortG;
addpath ./forked_MESHND;
addpath ./Morton;
addpath ./Hilbert;
addpath ./Plotting;
flag_save_files = 0;

q = 2;
n = 2^q;
h2=(n+1)*(n+1);
N = n

K1D = h2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);   % 1d Poisson matrix (negative Laplacian)

A = full(K1D);
b = 100*ones(N,1)
% x = K1D\b;

% Ab = [ -50  0;
% 	   0  -50];

% Bb = [  25  0;
% 	   25  25];

% Cb = [  25  25;
% 	    0  25];

% Db = [ -50  0;
% 	   0  -50];

% A11 = Db-Cb*inv(Ab)*Bb;
% A22 = Ab-Bb*inv(Db)*Cb;

% [100;100]-Cb*inv(Ab)*[100;100];
% [100;100]-Cb*inv(Ab)*[100;100];

P = [-31.25 0   	0  		0;
	   	0  	-20.833 0   	0;
	   	0  	0   	-20.833 0;
	   	0  	0		0  		-31.25];

b = [250;250;250;250]

inv(P)*b

eigA  = eig(A);
eigPA = eig(P*A);

condA  = cond(A);
condPA = cond(P*A)

% P*A;
% inv(A)*A;

% p = [1,3,2,4]; % red/black
% p = [2,4,1,3]; % black/red
% full(K1D(p,p));

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