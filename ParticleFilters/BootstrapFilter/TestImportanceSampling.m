x_min = 0;
x_max = 2;
proposal_distribution = Uniform(x_min, x_max);

mean = 1.0;
variance = 0.02;
target_distribution = Gaussian(mean, variance);
target_pdf = @(x) target_distribution.pdf(x);

n_realizations = 200;

% Draw n realizations of the proposal distribution and weigh them according to the target pdf
[x, weights] = ImportanceSampling_WeightedSamples(target_pdf, proposal_distribution, n_realizations)

% Convert these draw weighted samples into a Histogram
hist = Histogram(x_min, x_max, 25, x, weights);

figure(1);
subplot(2,1,1);
stem(x, n_realizations*weights/2);
hold on;
target_distribution.plot();
hold off;

subplot(2,1,2);
hist.plot();