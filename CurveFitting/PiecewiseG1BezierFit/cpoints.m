function C = cpoints(x)

[P,ang,dt] = ktangdt(x); % Separate vector x.
C = ctpts(P,ang,dt); % Call to compute the control points.