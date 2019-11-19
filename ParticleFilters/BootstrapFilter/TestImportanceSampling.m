proposal_distribution = Uniform(0, 2);

mean = 0.5;
variance = 0.02;
target_distribution = Gaussian(mean, variance);
target_pdf = @(x) target_distribution.pdf(x);

n_realizations = 200;

[x, weights] = ImportanceSampling(target_pdf, proposal_distribution, n_realizations)

fig = figure(2);
stem(x, n_realizations*weights/2);
hold on;
target_distribution.plot(fig);
hold off;