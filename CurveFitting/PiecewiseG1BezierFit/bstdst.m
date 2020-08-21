function bdt = bstdst(id,Q,P,ang,k)

% This function finds the optimum distances for control point
% placement along the segments of a curve. The applicable points
% from Q, the two knots, and two angles for each segment are
% passed to opdist.m through "fmins". It was written by E. J. Lane.

% opts = [0, 0.01, 0.01]; % Control parameters for "fmins".

options.Display = 'off';  % ECR, to match obselete fmins params above
options.TolX    = 0.01;
options.TolFun  = 0.01;

n = length(id); bdt=[];

for i = 1:n
   bdt(:,i) = fminsearch(@opdist, id(:,i), options,...
       Q(:,k(i):k(i+1)), P(:,i:i+1), ang(i:i+1));
end