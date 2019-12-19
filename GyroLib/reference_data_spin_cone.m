function [dWb, q_ref, dt, q_init, Fb, Mb, fn, mn, bw] = reference_data_spin_cone
% [dWb, q, dt, q_init, fb, mb, fn, mn, bw] = reference_data_spin_cone
% Creates reference data for AHRS algorithm simulation
% Reference movement of the body is defined by the SPIN-CONE truth model -
% body spinning at a fixed magnitude rotation rate and whose spin axis is
% rotating at a fixed precessional rate
% Ref: Paul G. Savage, Strapdown System Performance Analysis
%
%   Output arguments:
%   dWb -  Gyroscopes measurements (integral of the angular rate over 
%   the computer cycle)
%   q  - Reference attitude quaternion readings
%   dt - time interval between measurements (computer cycle)
%   q_init - Initial reference quaternion
%   fb  - Accelerometers measurements
%   mb  - Magnetometer measurements
%   fn - Gravity vector in navigation frame
%   mn - Magnetic field vector in navigation frame 
%   bw - Gyroscopes model biases

%%
clear;
clc;

%%
rng('shuffle');

%% Simulation parameters
%Simulation time, sec.
tsim = 360;
%Simulation time step, sec.
dt = 1e-2;

%% Kinematic parameters
%Perecessional rate, rad/sec
wc = 1e-1;
%Body spin rate, rad/sec
ws = 1e0;
%Angle between the precessional axis and spin axis - "cone angle"
beta = pi/4;

%% Sensors Noises
nw = 1e-8; %gyroscopes
nf = 1e-6; %accelerometers
nm = 1e-6; %magnetometers

%% Sensors biases
bdw = [0e0; 0e0; 0e0]; %gyroscopes
bfb = [0e0; 0e0; 0e0]; %accelerometers
bmb = [0e0; 0e0; 0e0]; %magnetometers

%% Biases random walk 
rwn = 1e-8; %gyroscopes

%% Navigation frame vectors 
fn = [0;0;-1]; %"gravity" in n-frame
mn = [1;0;0];  %"magnetic field" in n-frame

%%
dWb = zeros(tsim/dt,3);
q_ref   = zeros(tsim/dt,4);
angles = zeros(tsim/dt,3);
Fb = zeros(tsim/dt,3);
Mb = zeros(tsim/dt,3);
bw = zeros(tsim/dt,3);

%% Simulation
time = 0;
cnt = 0;
[iWb_, Cbn] = spin_cone( wc, ws, beta, time );
q_init = dcm_quat(Cbn);
while (time < tsim)
    time = time+dt;
    cnt = cnt+1;
    bdw = bdw+randn(3,1)*rwn;
    [ iWb, Cbn, rz, ry, rx ] = spin_cone( wc, ws, beta, time );
    da = iWb-iWb_;
    dWb(cnt,:) = da+randn(3,1)*sqrt(nw)+bdw;
    q_ref(cnt,:) = dcm_quat(Cbn);
    angles(cnt,:) = [rz, ry, rx];
    Fb(cnt,:) = Cbn'*fn+randn(3,1)*sqrt(nf)+bfb;
    Mb(cnt,:) = Cbn'*mn+randn(3,1)*sqrt(nm)+bmb;
    bw(cnt,:) = bdw;
    iWb_ = iWb;
end

