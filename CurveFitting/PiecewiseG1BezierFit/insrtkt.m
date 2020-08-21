function [xi,nk] = insrtkt(seg,h,x0,k,Q)

% function [xi,nk] = insrtkt(seg,h,x0,k). This function receives the
% segment number, seg, to have the knot inserted, the position along
% the segment, h, where it will be inserted, and the vector, x0, of
% parameters for the curve, and k, the knot positions. It inserts a
% new knot on the segment called for and then returns the new vector
% of parameters for the curve and knot positions. Note: the curve
% will remain the same, the polygon will be changed. This was
% written by E. J. Lane.

[P,ang,dt] = ktangdt(x0); % Separates xO into its subcomponents.

q = length(k);

% Computes the control points for the affected segment.
Cseg = ctpts(P(:,seg:seg+1),ang(seg:seg+1),dt(:,seg));

z=fndpts(Cseg,h); 	% Call to compute new control points for the
% segment where the knot is inserted.

xs=z(1,:); ys=z(2,:); 	% Separates the new segment's control
% points into their x and y components.

dx   = diff(xs); dy=diff(ys); % Finds the intercomponent differences.
angs = atan2(dy,dx); 	      % Computes the angles for the tangent vectors.

% Computes distances for control point locations.
dl = sqrt(dx(1)^2 + dy(1)^2);
de = sqrt(dx(6)^2 + dy(6)^2);
dm = sqrt(dx(4)^2 + dy(4)^2);
dn = sqrt(dx(3)^2 + dy(3)^2);

% Inserts new knot into knot component vector.
Pnew = [P(:,1:seg) z(:,4) P(:,seg+1:length(P))];

% Inserts new angles into tangent angles component vector.
angnew = [ang(1:seg) angs(4) ang(seg+1:length(ang))];

% Inserts new distances into distance component vector.
dtnew = [dt(:,1:seg-1) [dl dm;dn de] dt(:,seg+1:length(dt))];

dv = Q(:,k(seg):k(seg+1)) - z(:,4)*ones(1,k(seg+1)-k(seg)+1);
ds = dv.*dv; dq = sum(ds); 
[dmin,knew] = min(dq);

% With previous 2 lines finds the new knot's position.
ink = k(seg) + knew - 1;

% New knot sequence.
nk = [k(1:seg) ink k(seg+1:q)];

% Assembles the new components vector for the
% parameters for the curve
xi = [Pnew(1,:) Pnew(2,:) angnew dtnew(1,:) dtnew(2,:)];

