function dt = distEJL(P)

% function dt = dist(P). This function computes the initial
% distances from the knot points to their adjacent control
% points for the initial guess curve. It returns the vector
% of distances to iguess.m. The function was written by
% E. J. Lane.

t  = length(P);
d1 = P(:,1:t-1) - P(:,2:t); 	% Calculates inter-knot
				% x and y difference values.

d2 = sqrt(sum(d1.^2))/3; 	% Computes the initial distances.
dt = [d2;d2]; 			% Assembles the vector of distances.
