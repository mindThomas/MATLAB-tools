function se = objf2(x,Q,t,k)

% function se = objf2(x,Q,t,k). This is objective function that 
% will be minimized by "fmins”. The input arguments are the vector x of
% parameters for the curve, data points Q, a toggle t if a new knot has
% been inserted or one removed, and the knot sequence, k. The output is
% the sum from the function “sod” plus the the distance square from the
% first and last data points to the first and last knot points, respec0-
% tively. This was written by M. R. Holmes and revised by E. J. Lane.

global dpkpc

% LLoop to change dpkpc if a knot was inserted or removed.
if t == 1
	dpkpc = k; t = n;
	global dpkpc
end

if t == 0
	global dpkpc
	[r,s] = size(Q);
	[P,ang,dt] = ktangdt(x); % Call to separate x into its subcomponents.
	C = ctpts(P,ang,dt);     % Call to compute control points.
	dpkpc = newk(Q, P); 	 % Calls  function computes the
				 % new dividing point positions .
	m =  length(x);
	n = round(m/5);
	fp = P(:,1) - Q (:, 1);  % Compute the distance squared
	lp = P(:,n) - Q(:,s); 	 % from the first and last data points to
				 % the first and last knot points, respectively.

	% Calls the function that computes the sums of the square
	% of the distances from the data points to the nearest point
	% on the cubic segment.
	se = sod(C,Q,dpkpc) + fp'*fp + lp'*lp;
end
end
