function asset_NestedDissection_Prototype
clear all; close all; clc; format compact;

tic
q = 5; %size of the problem

plotDense_flag = 0;
plotSpy_flag = 1;
plotSol_flag = 0;

n = power(2,q);

% 2D and 3D meshes
stencils = [5 9 7 27];
s = 1; %2D
% s = 3; %3D

nn = [n n n n]; mm = [n n n n]; kk = [1 1 n n]; 
m = mm (s); n = nn (s); k = kk (s);

Gnew = reshape (1:(m*n*k), n, m, k)'; % List of natural order (row or colum-major)
A = -meshsparse (Gnew, stencils (s)); % FD matrix

% Ordering schemes:
p = nd2(Gnew); % Nested ordering
p0 = p-ones(size(p));
% save p0.mat p0 -ascii
size(p)

% p = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16] % Natural
% p = [1,9,13,5,3,4,11,15,12,16,7,8,2,6,10,14]; % Nested ordering
% p = [1 2 3 4 9 10 11 12 5 6 7 8 13 14 15 16]; % Red/Black ordering
% p = [7,13,5,4,9,16,3,6,12,11,15,2,10,8,14,1]; % Inverse METIS
% p = [16,12,7,4,3,8,1,14,5,13,10,9,2,15,11,6]; %  METIS
% C = A (p,p); % Permuted sparse matrix
C = A (p,p);
% full(C)

if plotSpy_flag ==1
    figure(1)
    subplot(1,2,1); do_spy (A); title ('natural ordering');
    subplot(1,2,2); do_spy (C); title ('permuted'); 
end

% return

% Comment out solve:
% % figure(2)
% % Set-up a test RHS
% b = ones(n*m*k,1);
% xOr = A\b; % Solve the system

% % Solving by LU decomposition
% [L,U] = lu(C);
% xPer = U\(L\b);

% % Get reversal of a sort (Permute)
% % http://blogs.mathworks.com/loren/2007/08/21/reversal-of-a-sort/#7
% inv_p(p) = 1:(m*n*k); % inv_p has the reverse indices
% sorted_xPer = xPer(inv_p);


% if plotSol_flag ==1
%     figure(2)
%     plot(xOr,'o','Color','blue');
% 	hold on;
%     plot(sorted_xPer,'*','Color','red'); title('Superposing solutions with different orderings');
% end

% % Show permutation vectors
% if plotDense_flag ==1
%     oriA = full(A)
%     perC = full(C)
%     p'-ones(size(p))'
%     inv_p;
% end
% size(p);

% toc
% % err = L*U*x-b;  % Numerical error of decomposition
% % err = xOr-sorted_xPer; % Numerical error against exact solution
% % printf('Solution error is: %1.1f ', norm(err));

