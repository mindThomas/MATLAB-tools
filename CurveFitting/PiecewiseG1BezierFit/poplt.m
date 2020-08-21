function pop = poplt(x,Q)


% function pop = poplt(x,Q). This function picks out subcomponents
% of the vector x of curve parameters. It calls the function that
% computes the control points. It then calls for a plot of the curve
% its polygon, and the data points. This was written by M. R. Holmes
% and revised by E. J. Lane.

[P,ang,dt] = ktangdt(x); % Separate vector x.

C = ctpts(P,ang,dt); 	 % Call to compute the control points.

pop = pltC(C,Q,P); 		 % Call to plot the curve, polygon, and data points.