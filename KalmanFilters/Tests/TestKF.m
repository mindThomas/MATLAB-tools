scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, '../'));
addpath(fullfile(scriptDir, '../MotionModels'));
addpath(fullfile(scriptDir, '../MeasurementModels'));

%% Create motion model
ts = 0.1;
[f, Fx, Fu, Fq] = CoordinatedTurnModel_Discrete(ts);

% x = [ x, y, v, phi, omega ]
x0 = [ 0.2, 0, 0.2, deg2rad(25), deg2rad(25) ]';
u0 = zeros(size(Fu(zeros(99,1),zeros(99,1),zeros(99,1)),2),1);
q0 = zeros(size(Fq(zeros(99,1),zeros(99,1),zeros(99,1)),2),1);

% Define process covariance
sigma_q_v = 0.1;
sigma_q_omega = 0.01; 
Q = diag([sigma_q_v^2, sigma_q_omega^2]);

%% Create measurement model
[h, Hx, Hr] = PositionSensor(1,2); % specify where in the state vector that (x,y) is located

% z = [ z_x, z_y ]
r0 = zeros(size(Hr(zeros(99,1),zeros(99,1)),2),1);

% Define measurement covariance
sigma_r_x = 1;
sigma_r_y = 1; 
R = diag([sigma_r_x^2, sigma_r_y^2]);

%% Initialize Kalman filter
kf = KF;

% Define initial process covariance
% x = [ x, y, v, phi, omega ]
sigma_xy = 1;
sigma_v = 0.1;
sigma_phi = 0.01;
sigma_omega = 0.01;
P0 = diag([sigma_xy^2, sigma_xy^2, sigma_v^2, sigma_phi^2, sigma_omega^2]);

x_init = x0;
x0 = zeros(size(x0));
kf = kf.init_discrete_function_jacobians(...
                                f, Fx, Fu, Fq, x0, u0, q0, Q, ...  % process model
                                h, Hx, Hr, R, r0, ...              % measurement model
                                P0);
kf.x = x_init;                            
               
time = 0;                    
pos = kf.x(1:2)';
x = x_init;
for (i = 1:100)    
    time(end+1,1) = time(end,1) + ts;
    %kf = kf.predict();
    x = f(x,u0,q0);
    pos(end+1,:) = x(1:2)';
end

figure(1);
%plot(time, pos(:,1));
plot(pos(:,1), pos(:,2));
axis equal;