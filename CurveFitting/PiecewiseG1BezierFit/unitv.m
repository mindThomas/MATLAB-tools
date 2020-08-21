function uv = unitv(Q,k)

% function uv = unitv(Q,k). This function takes data points,
% Q, and the position of knotpoints, k, as input variables.
% It uses chord length parameterization to fit a parametric
% quadratic curve to five data points. The unit tangent vec-
% tors are approximated by the unit tangent vectors for these
% quadratic functions. It returns the set of unit tangent
% vectors in the direction of the knot points. It was written
% by M. R. Holmes

[r, m] = size (Q);
n  = length(k);
for j = 1:n 	% Loop to index knot positions.
    if j == 1
        k(j) = 1; kt = 1;
    elseif j == n
        k(j) = m-4; kt = 5;
    else
        k(j) = k(j)-2; kt = 3;
    end
    
    % Extracting the knot point and four adjacent points.
    x = Q(1, k(j):k(j)+4)';
    y = Q(2, k(j):k(j)+4)';
    
    xd = diff(x); yd= diff(y); % Get chord length.
    d  = sqrt(xd.*xd + yd.*yd);
    t(1) = 0;
    t(2) = d(1);
    t(3) = t(2) + d(2);
    t(4) = t(3) + d(3);
    t(5) = t(4) + d(4);
    
    c = [ones(5,1) t' (t.*t)'] \ [x y];
    u = c(2, :) + 2*c(3,:)*t(kt);
    u = u/norm(u);
    uv(:,j) = u'; 	% Approximation of unit tangents
                    % by unit tangent to quadratic.
end