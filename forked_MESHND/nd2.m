function p = nd2 (G)
%ND2 p = nd2 (G) permutes a 2D or 3D mesh G.
% Compare with nestdiss which uses p as a scalar offset and returns a modified
% mesh G that corresponds to Gnew in meshnd.  Here, the scalar offset p in
% nestdiss is not needed.  Instead, p is a permutation, and the modified mesh
% Gnew is not returned.

[m n] = size (G);
% printf('Empezo: (m,n) = (%d,%d)',m,n);
G;

if ( m <= 2) % G is small; do not cut it
    % disp('caso1');
    p = G (:);

% elseif n >= m % cut G along the middle column, cutting n in half
elseif n >= m % cut G along the middle column, cutting n in half
    % disp('caso3');
    s = ceil(n/2);
    c3_p1 = G(:,1:s-1);
    c3_p2 = G(:,s+1:n);
    c3_middle = G(:,s);
    p = [ (nd2(c3_p1))  ; (nd2(c3_p2)) ; c3_middle(:) ];

else % cut G along the middle row, cutting m in half   
    % disp('caso4');
    s = ceil(m/2);
    c4_p1 = G(1:s-1,:);
    c4_p2 = G(s+1:m,:);
    c4_middle = G(s,:);
    p = [(nd2(c4_p1)) ; (nd2(c4_p2)) ; c4_middle(:)];

end