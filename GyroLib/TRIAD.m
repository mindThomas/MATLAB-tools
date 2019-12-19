function Cbn = TRIAD(fb, mb, fn, mn)
% Cbn = TRIAD(fb, mb, fn, mn)
% Function implements  TRIAD algorithm using measurements
% from three-component accelerometer with orthogonal axes and vector
% magnetometer
%
%   Input arguments:
%   fb  - Acceleration vector in body frame [3x1]
%   mb  - Magnetic field vector in body frame [3x1]
%   fn  - Gravity vector in navigation frame [3x1]
%   mn  - Magetic field vector in navigation frame [3x1]
%
%   Output arguments:
%   Cbn - estimated Direction Cosines Matrix (DCM)

W1 = fb/norm(fb);
W2 = mb/norm(mb);

V1 = fn;
V2 = mn;

Ou1 = W1;
Ou2 = cross(W1,W2)/norm(cross(W1,W2));
Ou3 = cross(W1,cross(W1,W2))/norm(cross(W1,W2));

R1 = V1;
R2 = cross(V1,V2)/norm(cross(V1,V2));
R3 = cross(V1,cross(V1,V2))/norm(cross(V1,V2));

Mou = [Ou1, Ou2, Ou3];
Mr = [R1, R2, R3];

%TRIAD DCM
Cbn = Mr*Mou';

end