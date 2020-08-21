function error = err(x0,Q,k)

% function error = err(xO,Q,k). This function takes a composite
% vector of curve parameters x0, separates them and computes the
% control points for the curve. It then computes the sum of the
% error between the curve and the data points in Q. It was
% written by E. J. Lane.

[P,ang,dt] = ktangdt(x0); % Call to separates x0
			% into its subcomponents.

C = ctpts(P,ang,dt); 	% Call to compute control points.

error = sod(C, Q, k); 	% Call to compute distance
			% error for the curve.
