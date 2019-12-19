function [q, P, bw] = ahrs_quat(q, P, bw, dwb, fb, mb, fn, mn, dt)
% [q, P, bw] = ahrs_quat(q, P, bw, dwb, fb, mb, fn, mn, dt)
% Implements the quaternion AHRS algorithm usin measurements
% from the three-axis gyroscope and accelerometer,
%  vector magnetometer
%
%   Input arguments:
%   q   - Attitude quaternion [1x4]
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
%   q  - Updated attitude quaternion [1x4]
%   P  - Updated Kalman filter covarince matrix [6x6]
%   bw - Updated gyroscopes biases vector [3x1]

%% Correct attitide increments  with estimated biases
dwb  = dwb-bw*dt;

%% Attitude Quaternion Mechanization
%Quaternion from k+1 body axes to k body axes
gamma1 = dwb(1);
gamma2 = dwb(2);
gamma3 = dwb(3);
gamma = norm(dwb);
lambda0 =cos(gamma/2);
lambda1 = gamma1*sin(gamma/2)/gamma;
lambda2 = gamma2*sin(gamma/2)/gamma;
lambda3 = gamma3*sin(gamma/2)/gamma;
if (gamma~=0)
    lambda = [lambda0, lambda1, lambda2, lambda3];
else
    lambda = [1, 0, 0, 0];
end
%Update q for body motion
q = quat_mult(quat_conj(lambda),q);

%Normalize quaternion
q = q/sqrt(q*q');

%% Kalman predict
%System dynamics matrix F and system noise covariance matrix Q
%Continuous-time system matrix
A = zeros(6,6);
A(1:3,4:6) = quat_dcm(q);
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
%Estimated measurements:
%"Gravity" vector estimate in navigation frame
fn_hat = quat_rot(q,fb);
%"Magnetic field" vector estimate in navigation frame
mn_hat = quat_rot(q,mb);

%Measurement vector
v = zeros(6,1);
v(1:3,1) = fn-fn_hat;
v(4:6,1) = mn-mn_hat;

%Measurement matrix
H = zeros(6,6);
H(1:3,1:3) = skew(fn_hat);
H(4:6,1:3) = skew(mn_hat);

%Measurement noise covariance matrix
R = diag([1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1]);

%% Kalman Update
I = eye(6);
S = H*P*H+R;
K = (P*H')/S;
P = (I-K*H)*P*(I-K*H)'+K*R*K';
x = K*v;

%% Correct Attitude Quaternion
%Error quaternion
qe = [1, x(1)/2, x(2)/2, x(3)/2];
qe = qe/sqrt(qe*qe');
%Correct esimated attitude quaternion
q = quat_mult(q, qe);

%% Update gyro bias estimate
bw = bw+x(4:6,1);

end


