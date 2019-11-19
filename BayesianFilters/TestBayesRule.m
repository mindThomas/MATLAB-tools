%% Filtering equations
% Prior distribution   (prediction step)
% p(x[k] | y[1:k-1]) = integral p(x[k] | x[k-1]) * p(x[k-1] | y[1:k-1]) dx[k-1]
%   where p(x[k-1] | y[1:k-1]) is the previous posterior distribution
%
% Posterior distribution   (update step)
% p(x[k] | y[1:k]) = p(y[k] | x[k]) * p(x[k] | y[1:k-1])   /   p(y[k])
% Which is proportional to:
% p(x[k] | y[1:k]) ~ p(y[k] | x[k]) * p(x[k] | y[1:k-1])
%   where p(y[k] | x[k]) is the likelihood coming from the measurement model

%% Test setup
% Lets start with a Gaussian centered at 0 as our previous posterior distribution
xmin = -6; xmax = 6;
x_true = -4;
mu = -4;
sigma = 1;
initial = Gaussian(mu, sigma^2);
posterior = initial;

%%
figure(1);
subplot(3,2,1);
posterior.plot();
title('Previous posterior');
xlim([xmin, xmax]);

% Lets define propagation model as shifting the mean and increasing the
% variance according to:
% x[k] = x[k-1] + 1.0 + q
%   where q ~ N(0, 0.5)
% This model thus defines:
%   p(x[k] | x[k-1])
propagation = ConditionalGaussian(1.0, 1, 0.2);
x_true = x_true + 1.0;
subplot(3,2,2);
propagation.plot();
title('Propagation distribution');
xlim([xmin, xmax]);
ylim([0, 0.6]);

% Prediction step
% First we form a joint distribution as:
% p(x[k], x[k-1] | y[1:k-1]) = p(x[k] | x[k-1]) * p(x[k-1] | y[1:k-1])
joint = propagation.join(posterior);
% Next we marginalize out the previous state, x[k-1], to get:
% p(x[k] | y[1:k-1]) = integral p(x[k], x[k-1] | y[1:k-1]) dx[k-1]
prior = joint.marginalize(2);
subplot(3,2,3);
prior.plot();
title('Prior distribution');
%xlim([xmin, xmax]);
ylim([0, 0.6]);

% ### Bayes rule is applied in the Update step ###
% Lets first define our measurement model as directly measurement the state
% plus some noise
% z[k] = x[k] + r
%   where r ~ N(0, 1.0)
% This model thus defines:
%   p(y[k] | x[k])
var_meas = 0.5;
measurement = ConditionalGaussian(0.0, 1, var_meas);
subplot(3,2,4);
measurement.plot();
title('Measurement distribution');
xlim([xmin, xmax]);
ylim([0, 1.4]);

% Update step
% First we form a joint distribution as:
% p(y[k], x[k] | y[1:k-1]) = p(y[k] | x[k]) * p(x[k] | y[1:k-1])
joint = measurement.join(prior);
% Which can be split using the product rule of Conditional probability
% p(y[k], x[k] | y[1:k-1]) = p(x[k] | y[1:k-1], y[k]) * p(y[k])
% Assuming p(y[k]) to be uniform, we can see this as a normalization factor
% Thus we can get the final posterior by splitting the joint distribution
% into its' conditional counterpart and evaluating it at the measured value
y = x_true + sqrt(var_meas)*randn(1,1); % simulate a realization of a measurement
posterior = joint.conditional(1, y)
subplot(3,2,5);
posterior.plot();
hold on;
vline(y, 'g', 'Measurement');
vline(x_true, 'r', 'True state');
hold off;
title('Posterior distribution');
%xlim([xmin, xmax]);
ylim([0, 0.7]);