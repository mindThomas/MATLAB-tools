function GOC = globop(xi,Q,t,k)

% function GOC = globop(xi,Q,t,k) This function returns a vector
% GOC of parameters: knot points, P, angles, ang, and distances,
% dt, for a globally optimized Bezier curve. Its inputs are
% the curve parameters in vector xi , data points, Q, toggle, t,
% "l" if a knot was inserted or removed, "0" otherwise , and the
% knot sequence, k. The MATLAB routine "fmins" optimizes function
% objf2.m which computes the sum of the distances between the data
% points and their closest point on the curve. It was written by
% E. J. Lane.

% GOC = fmins('objf2', xi, (0, 0.01, 0.01], [], Q, t, k);

options.Display = 'off';  % ECR, to match obselete fmins params above
options.TolX    = 0.01;
options.TolFun  = 0.01;

GOC = fminsearch(@objf2, xi, options, Q, t, k);