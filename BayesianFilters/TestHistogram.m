mu = 0;
cov = 1;
dist = Gaussian(mu, cov);
x_min = -6;
x_max = 6;
num_bins = 8*(x_max-x_min) + 1;
hist = Histogram(x_min, x_max, num_bins, @(x)(dist.pdf(x)));
hist.plot();

%%
f = @(x) x.^2;
hist2 = hist.propagate(f);
hold on;
hist2.plot();
hold off;

%%
mu = [0;0];
cov = diag([1;2]);
dist = Gaussian(mu, cov);

x_min = [-2*3.5; -2*3.5];
x_max = [2*3.5; 2*3.5];
num_bins = [4*7; 4*7];
hist = Histogram(x_min, x_max, num_bins, @(x)dist.pdf(x));
hist.plot();

