function [Cbn, P, bw] = ahrs_dcm(Cbn, P, bw, dwb, fb, mb, fn, mn, dt)
% [Cbn, P, bw] = ahrs_dcm(Cbn, P, bw, dwb, fb, mb, fn, mn, dt)
% Implements the direction cosine matrix AHRS algorithm using 
% measurements from the three-axis gyroscope and accelerometer,
%  vector magnetometer
%
%   Input arguments:
%   Cbn - Direction cosine matrix [3x3]
%   P   - Kalman filter covarince matrix [6x6]
%   bw  - Gyroscopes biases vector [3x1]
%   dwb - Integral of angular rate vector over the computer cycle [3x1]
%   fb  - Acceleration vector in body frame [3x1]
%   mb  - Magnetic field vector in body frame [3x1]
%   fn  - Gravity vector in navigation frame [3x1]
%   mn  - Magetic field vector in navigation frame [3x1]
%   dt  - Computer cycle, sec.
%
%   Output arguments:
%   Cbn  - Updated direction cosine matrix [3x3]
%   P    - Updated Kalman filter covarince matrix [6x6]
%   bw   - Updated gyroscopes biases vector [3x1]


%% Correct attitide increments  with estimated biases
dwb = dwb-bw*dt;

%% Attitude DCM Mechanization
%Cbb DCM from k+1 body axes to k body axes
rot_norm=norm(dwb);
sr_a=1-(rot_norm^2/6)+(rot_norm^4/120);
sr_b=(1/2)-(rot_norm^2/24)+(rot_norm^4/720);
Cbb=eye(3)+sr_a*skew(dwb)+sr_b*skew(dwb)*skew(dwb);
%Update Cbn for body motion
Cbn=Cbn*Cbb;

%% Kalman predict
%System dynamics matrix F and system noise covariance matrix Q
%Continuous-time system matrix
A = zeros(6,6);
A(1:3,4:6) = -Cbn;
%State transition matrix
F = eye(6)+A*dt+A*A*dt*dt/2;
%Noise-input mapping matrix
G = zeros(6,6);
G(1:3,1:3) = eye(3);
G(4:6,4:6) = eye(3);
%Gyro errors noise
ng = 1e-5;
%Gyro bias noise
ngb = 1e-16;
%System noise
Qn = diag([ng, ng, ng, ngb, ngb, ngb]);
%Trapezioidal integration
Q = 1/2*(F*G*Qn*G'+G*Qn*G'*F')*dt;

%Covariance predict
P = F*P*F'+Q;

%% Measurements
%Estimated measurements
%"Gravity" vector estimate in navigation frame
fn_hat = Cbn*fb;
%"Magnetic field" vector estimate in navigation frame
mn_hat = Cbn*mb;

%Measurement vector
v = zeros(6,1);
v(1:3,1) = fn_hat-fn;
v(4:6,1) = mn_hat-mn;

%Measurement matrix
H = zeros(6,6);
H(1:3,1:3) = skew(fn_hat);
H(4:6,1:3) = skew(mn_hat);
R = diag([1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1]);

%% Kalman Update
I = eye(6);
S = H*P*H'+ R;
K = (P*H')/S;
P = (I-K*H)*P*(I-K*H)' + K*R*K';
x = K*v;

%% Correct Attitude DCM
E = eye(3)+skew(x(1:3,1));
Cbn = E*Cbn;
%Normalize Cbn matrix
delta_12 = Cbn(1,:)*Cbn(2,:)';
Cbn(1,:) = Cbn(1,:) - 1/2*delta_12*Cbn(2,:);
Cbn(2,:) = Cbn(2,:) - 1/2*delta_12*Cbn(1,:);
Cbn(3,:) = cross(Cbn(1,:),Cbn(2,:));
Cbn(1,:) = Cbn(1,:)./norm(Cbn(1,:));
Cbn(2,:) = Cbn(2,:)./norm(Cbn(2,:));
Cbn(3,:) = Cbn(3,:)./norm(Cbn(3,:));

%% Update gyro bias estimate
bw = bw+x(4:6,1);

end
