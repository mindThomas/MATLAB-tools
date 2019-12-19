% Bearing only tracking problem with constant velocity motion model in 2D

T = 0.1; % sample time
% Linear states: 2D Velocity
% x1[k] = x1[k-1] + qx1[k-1]
% x2[k] = x2[k-1] + qx2[k-1]
% Motion model has no dependence on non-linear states
f = @(y) zeros(2,1);
Fx = @(y) eye(2);
velocity_cov = 1 * diag([1; 1]);

velocity_to_position = [T*eye(2); T^2/2*eye(2)];
cov_xy = velocity_to_position * velocity_cov * velocity_to_position';

% Non-linear states: 2D Position (since they enter the measurement model non-linearly)
% y1[k] = y1[k-1] + T*x1[k-1] + qy1[k-1]
% y2[k] = y2[k-1] + T*x2[k-1] + qy2[k-1]
% Motion model depends on both non-linear and linear states
g = @(y) y;
Gx = @(y) T*eye(2);

% Measurement model: Bearing
% z[k] = atan2(y2[k], y1[k]) + r[k]
h_bearing = @(y) atan2(y(2), y(1));
h_range = @(y) sqrt(y(1)^2+y(2)^2);
Hx_bearing = @(y) zeros(1,2);
Hx_range = @(y) zeros(1,2);
h = @(y) [h_bearing(y); h_range(y)];
Hx = @(y)[Hx_bearing(y); Hx_range(y)];
cov_bearing = deg2rad(1)^2; % 1 deg standard deviation
cov_range = 0.1^2; % 3 meter standard deviation
% Measurement model conditional PDF
meas_bearing_pdf = @(z, y) 1/sqrt((2*pi)*det(cov_bearing)) * exp(-1/2 * (z-h_bearing(y))' * inv(cov_bearing) * (z-h_bearing(y)));
meas_range_pdf = @(z, y) 1/sqrt((2*pi)*det(cov_range)) * exp(-1/2 * (z-h_range(y))' * inv(cov_range) * (z-h_range(y)));
meas_pdf = @(z, y) meas_bearing_pdf(z(1), y) * meas_range_pdf(z(2), y);
cov_r = [cov_bearing, 0; 
         cov_range, 0];

%%
RBPF = RaoBlackwellizedPF(2, f, Fx, 2, g, Gx, cov_xy, h, Hx, cov_r, meas_pdf);

n_particles = 200;
initial_nonlinear_distribution = Uniform([0;0], [10;10]);
initial_linear_mean = [0;0];
initial_linear_covariance = 1 * eye(2);
RBPF = RBPF.init(n_particles, initial_nonlinear_distribution, initial_linear_mean, initial_linear_covariance)

figure(1);
[states, weights] = RBPF.getStatesAndWeights();
stem3(states(3,:), states(4,:), weights);
xlim([0,10]);
ylim([0,10]);

figure(2);
clf;
gaussians = RBPF.getLinearGaussians();
subplot(2,1,1);    
gaussians{1}.marginalize(1).plot()
subplot(2,1,2);    
gaussians{1}.marginalize(2).plot()
for (i = 2:length(gaussians))
    subplot(2,1,1);
    hold on;
    gaussians{i}.marginalize(1).plot()
    hold off;
    subplot(2,1,2);
    hold on;
    gaussians{i}.marginalize(2).plot()
    hold off;    
end

true_state = [1; -0.5;
              0.1; 5];

%%
%for (i = 1:100)
true_state(3:4) = true_state(3:4) + T*true_state(1:2);
measurement_bearing = atan2(true_state(4), true_state(3));
measurement_range = sqrt(true_state(3) + true_state(4)^2);
measurement = [measurement_bearing; measurement_range];
%measurement = deg2rad(90);
RBPF = RBPF.filter(measurement);

figure(1);
[states, weights] = RBPF.getStatesAndWeights();
stem3(states(3,:), states(4,:), weights);
xlim([0,25]);
ylim([-15,10]);
xlabel('x');
ylabel('y');

hold on;
plot(true_state(3), true_state(4), 'rv');
hold off;


%
figure(2);
clf;
gaussians = RBPF.getLinearGaussians();
subplot(2,1,1);    
gaussians{1}.marginalize(2).plot()
subplot(2,1,2);    
gaussians{1}.marginalize(1).plot()
for (i = 2:length(gaussians))
    subplot(2,1,1);
    hold on;
    gaussians{i}.marginalize(2).plot()
    hold off;
    subplot(2,1,2);
    hold on;
    gaussians{i}.marginalize(1).plot()
    hold off;    
end
%end

%%
[X,Y] = meshgrid(0:0.05:10, 0:0.05:10);
likelihood = [];
for (x = 1:size(X, 2))
    for (y = 1:size(Y, 2))
        pos = [x;y];
        likelihood(y,x) = meas_pdf(measurement, pos);
    end
end

surf(X,Y,likelihood);