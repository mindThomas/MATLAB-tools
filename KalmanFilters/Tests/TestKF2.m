scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, '../'));
addpath(fullfile(scriptDir, '../MotionModels'));
addpath(fullfile(scriptDir, '../MeasurementModels'));

%% Create motion model
ts = 0.1;
[f, Fx, Fu, Fq] = ConstantVelocity_Continuous();

% x = [ p, v ]
x0 = [ 0.2, 0.5 ]';
u0 = zeros(size(Fu(zeros(99,1),zeros(99,1),zeros(99,1)),2),1);
q0 = zeros(size(Fq(zeros(99,1),zeros(99,1),zeros(99,1)),2),1);

% Define process covariance
sigma_q_v = 0.05;
Q = sigma_q_v^2;

Qsim = 0.0001 * Q;
model = ContinuousMotionModelHandler(f, x0, Qsim);

%% Create measurement model
[h, Hx, Hr] = PositionSensor(1); % specify where in the state vector that (x,y) is located

% z = [ z_x, z_y ]
r0 = zeros(size(Hr(zeros(99,1),zeros(99,1)),2),1);

% Define measurement covariance
sigma_r_p = 0.01;
R = sigma_r_p^2;

Rsim = 1 * R;
meas = MeasurementModelHandler(h, Rsim);

%% Initialize Kalman filter
kf = KF;

% Define initial process covariance
% x = [ p, v ]
sigma_p = 0.1;
sigma_v = 0.01;
P0 = diag([sigma_p^2, sigma_v^2]);

x_init = x0;
x0 = zeros(size(x0));
% kf = kf.init_discrete_function_jacobians(...
%                                 f, Fx, Fu, Fq, x0, u0, q0, Q, ...  % process model
%                                 h, Hx, Hr, R, r0, ...              % measurement model
%                                 P0);
kf = kf.init_continuous_function_jacobians(...
                                f, Fx, Fu, Fq, x0, u0, q0, ts, Q, ...  % process model
                                h, Hx, Hr, R, r0, ...              % measurement model
                                P0);                       
kf.x = x_init;                            
               
time = 0;                    
true = x_init';
pred = kf.x';
est = kf.x';
measurements = x_init(1:2)';
for (i = 0:100)    
    time(end+1,1) = time(end,1) + ts;
    
    % Propagate model and store true position
    if (i == 50)
        model.x(2) = model.x(2) / 10;
    end
    model = model.step(ts);        
    true(end+1,:) = model.x';
    
    % Generate measurement
    z = meas.get(model.x);
    measurements(end+1,:) = z';
    
    % Kalman filter
    kf = kf.predict();    
    pred(end+1,:) = kf.x';
    %if (mod(i, 10) == 0)
        kf = kf.update(z);
    %end
    est(end+1,:) = kf.x';        
end

%%
figure(1);
subplot(2,1,1);
plot(time, true(:,1), time, pred(:,1), time, est(:,1));
hold on;
plot(time, measurements(:,1), 'r.');
hold off;
title('p'); legend('True', 'Prediction', 'Corrected', 'Measurements');
subplot(2,1,2);
plot(time, true(:,2), time, pred(:,2), time, est(:,2));
title('v'); legend('True', 'Prediction', 'Corrected', 'Measurements');