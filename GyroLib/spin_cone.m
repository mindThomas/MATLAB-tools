function [ iWb, Cbn, psi, theta, phi ] = spin_cone( wc, ws, beta, t )

%Attitude angles
phi = (ws-wc*cos(beta))*t;
theta = pi/2-beta;
psi = -wc*t;

%Direction cosine matrix
Cbn = zeros(3);
Cbn(1,1) = cos(theta)*cos(psi);
Cbn(1,2) = -cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi);
Cbn(1,3) = sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi);
Cbn(2,1) = cos(theta)*sin(psi);
Cbn(2,2) = cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi);
Cbn(2,3) = -sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi);
Cbn(3,1) = -sin(theta);
Cbn(3,2) = sin(phi)*cos(theta);
Cbn(3,3) = cos(phi)*cos(theta);

%Integrated angular rate
iWb = zeros(3,1);
if norm(ws) > 0
    iWb(1,1) = ws*t;
    iWb(2,1) = ((wc*sin(beta))/(ws-wc*cos(beta)))*cos(phi);
    iWb(3,1) = -((wc*sin(beta))/(ws-wc*cos(beta)))*sin(phi);
end

end

