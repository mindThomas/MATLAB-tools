function x = fndpts (z, h)

% function x = fndpts(z,h). The inputs are a vector z of control
% points for a segment of a curve and a step size h. The function
% separates the control points into their x and y components and
% then uses a de Casteljau or Chaikin scheme to compute new control
% points which will produce the same curve. It was written by
% E. J. Lane.

[m n]=size(z); M=zeros(n); N=zeros(n);
M(:,1) = z(l,:)'; 	% Separates the control points
N(:,1) = z(2,:)'; 	% x and y values.

% Loop which performs the computation
% of new control point x and y values.
for j = 2:n
    for i = j:n
        M(i, j) = M(i-1,j-1) + ((M(i,j-1) - M(i-1,j-1))*h);
        N(i, j) = N(i-1,j-1) + ((N(i,j-1) - N(i-1,j-1))*h);
    end
end

% Assembles the vector of new control points.
x = [diag(M)' (rot90(M(n,1:n-1)))'; diag(N)' (rot90 (N(n,1:n-1)))'];
