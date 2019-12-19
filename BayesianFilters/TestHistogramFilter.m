%% Test setup
% Lets start with a Gaussian centered at 0 as our previous posterior distribution
n = 200;
xmin = -8; xmax = 8;
mu = -4;
sigma = 1;
initial = Gaussian(mu, sigma^2);

% Lets define propagation model as shifting the mean and increasing the
% variance according to:
% x[k] = x[k-1] + 1.0 + q
%   where q ~ N(0, 0.5)
% This model thus defines:
%   p(x[k] | x[k-1])
propagation = ConditionalGaussian(1.0, 1, 0.2);

% Secondly we define our measurement model as directly measurement the state
% plus some noise
% z[k] = x[k] + r
%   where r ~ N(0, 1.0)
% This model thus defines:
%   p(y[k] | x[k])
var_meas = 0.5;
measurement = ConditionalGaussian(0.0, 1, var_meas);

% Run the Histogram filter
posterior = initial;

figure(1);
clf();
subplot(4,2,1);
posterior.plot();
xlim([xmin, xmax]); ylim([0 0.4]);
title('Analytically Propagation');
subplot(4,2,3);
propagation.plot();
hold on;
measurement.plot();
hold off;
legend('Propagation model', 'Measurement model');
xlim([xmin, xmax]); ylim([0 1]);

% Propagate analytically (since we are using Gaussians)
x_true = mu;
for (i = 1:10)
    % Visualize 10 predictions
    prior = propagation.join(posterior).marginalize(2);
    x_true = x_true + 1.0;
    
    subplot(4,2,5);
    hold on;
    prior.plot();    
    hold off;
    xlim([xmin, xmax]); ylim([0 1]);
    
    % And updates
    joint = measurement.join(prior);
    y = x_true + 0*sqrt(var_meas)*randn(1,1); % simulate a realization of a measurement
    posterior = joint.conditional(1, y);    
    %posterior = prior;
    
    subplot(4,2,7);
    hold on;
    posterior.plot();    
    hold off;
    xlim([xmin, xmax]); ylim([0 1]);
end
hold off;

%% Now try again with the Discrete Bayes filter using a Histogram/home/thomas/MATLAB/BayesianFilters/Distributions
hf = HistogramFilter(xmin, xmax, n, @(x,y)propagation.pdf(x,y), @(x,z)measurement.pdf(x,z), @(x)initial.pdf(x));

% Plot initial distribution
subplot(4,2,2);
hf.plot();
xlim([xmin, xmax]); ylim([0 0.4]*((xmax-xmin)/n));
title('Histogram (Discrete Bayes) Propagation');

% Plot propagation and measurement model
subplot(4,2,4);
Histogram(xmin, xmax, n, @(x)propagation.pdf(x,0)).plot();
hold on;
Histogram(xmin, xmax, n, @(x)measurement.pdf(x,0)).plot();
hold off;
legend('Propagation model', 'Measurement model');
xlim([xmin, xmax]); ylim([0 1]*((xmax-xmin)/n));

x_true = mu;
for (i = 1:10)
    % Visualize 10 predictions
    hf = hf.predict();
    x_true = x_true + 1.0;
    
    subplot(4,2,6);
    hold on;
    hf.plot();
    hold off;
    xlim([xmin, xmax]); ylim([0 1]*((xmax-xmin)/n));
    
    % And updates
    y = x_true + 0*sqrt(var_meas)*randn(1,1); % simulate a realization of a measurement
    hf = hf.update(y);
    
    subplot(4,2,8);
    hold on;
    hf.plot();
    hold off;
    xlim([xmin, xmax]); ylim([0 1]*((xmax-xmin)/n));
end