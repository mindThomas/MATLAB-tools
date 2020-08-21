function [xi,nk] = rmvkt(kt,x,k)

% function [xi,nk] = rmvkt(kt,x,k). This function takes inputs of
% which knot, kt, to remove, vector, x, of curve parameters, and kt
% sequence, k. It removes the knot, its angles, and its distances
% from the subcomponents of x, removes the index of the removed
% knot from k, then finds the knots that were adjacent to the one
% being removed, and constructs a new polygon for the "blended"
% curve segment. This was written by E. J. Lane.

[P,ang,dt] = ktangdt(x); % Separates the components of x.

n = length(P);

Pnew=[P(:,1:kt-1) P(:,kt+1:n)]; % Removes the knot.

m = length(ang);
angnew = [ang(:,1:kt-1) ang(:,kt+1:m)]; % Removes the knot's angles.

p = length(dt);
q = length (k);
nk = [k(1:kt-1) k(kt+1:q)]; % Get rid of removed knot in sequence.

% Computes the control points for the blended segment.
Cseg = ctpts(P(:,kt-1:kt+1),ang(kt-1:kt+1),dt(:,kt-1:kt));

xs = Cseg(1,:) ; ys=Cseg(2,:);	% Separates the x and y components.
dx = diff(xs); dy=diff(ys); 	% Gets the differences in the x's, y's.
dm = sqrt(dx(3)^2 + dy(3)^2);
dn = sqrt(dx(4)^2 + dy(4)^2); 	% Computes distances for the control
dmn = dm+dn; % points on the blended segment.

dt(1,kt-1) = dt(1,kt-1)*(dmn/dm);
dt(2,kt)   = dt(2,kt)  *(dmn/dn);

% Assembles the distances.

d1 = dt(1,:); d2= dt(2,:);
d11 = [d1(1:kt-1) d1(kt+1:p)];
d22 = [d2(1:kt-2) d2(kt:p)];
dtnew=[d11 ; d22];

% Assembles the composite vector of parameters for the curve.
xi = [Pnew(1,:) Pnew(2,:) angnew dtnew(1,:) dtnew(2,:)];
