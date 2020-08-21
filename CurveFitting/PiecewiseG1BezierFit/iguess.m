function [IG, k] = iguess(Q)

% function [IG,k]=iguess(Q). This routine takes a set of data points, Q;
% picks out subset of the data point for knot points, P; computes the
% position of the knot points,k; computes the initial distances, dt, to
% place the interior control points, C, which are also computed; compute
% the angles, ang, of the unit tangent vectors at each knot point; and
% assembles the vector IG of paramaters P, ang, and dt, for the curve.
% The routine returns the "vector" of parmeters and plots the curve,
% its polygon, and the data points in Q. It was written by M. R. Holmes
% and revised by E. J. Lane.

global dpkpc;
[r ,m] = size(Q);

disp('Give the number of knotpoints.');
n =input(' ');

disp('Type "1" for default knot position or "2" to input your own.');
h = input (' ');

if h == 1
    k = defk(m,n); % Calls for default knot position.
elseif h == 2
    disp('Input initial knot sequence as follows "[1 4 8 ... n]".')
    k =input(' ');
elseif (h ~= 1) || (h ~= 2)
    disp('Error! Start over and choose "1" or "2".') ,pause(2)
    iguess
end

dpkpc = k; 		% Position of knot points passed globally.
P   = knots(Q,k); 	% call to compute the knotpoints.
dt  = distEJL(P); 		% Call to compute the distance between
% successive ve knot points.

ang = tang(Q,k); 	% Call to compute the angles for
% the unit tangent vectors.

C   = ctpts(P,ang,dt);	% Call to compute the control points;
% for the curve.

pltC(C,Q,P); 		% Call to plot the initial guess curve,
% its control polygon, and points in Q.

% Assemble the composite vector of the initial guess curve parameters
IG = [P(1,:) P(2,:) ang dt(1,:) dt(2,:)];

