function A = morton(n)
% To create a Morton Scan order matrix
if n == 2
   A = [1 2; 3 4];
else
   B = morton(n/2);
   A = [B B+(n/2)^2; B+(n/2)^2*2 B+(n/2)^2*3];
end