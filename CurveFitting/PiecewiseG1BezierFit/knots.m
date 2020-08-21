function P = knots(Q,k)

% function P = knots(Q,k). This function takes data points Q
% and knot sequence vector k and picks out the knot points
% of the curve.

P= []; 
P= [P Q(:,k)];