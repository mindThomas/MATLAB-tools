function ang = tang(Q,k)

% function ang = tang(Q,k). This function computes the angles
% for the unit tangent vectors at the knot points.

u = unitv(Q,k); 		% Call to compute unit tangents.
ang = atan2(u(2,:), u(1,:));	% Converts tangents to angles.