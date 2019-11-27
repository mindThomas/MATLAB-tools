% Setup Non-linear filter benchmark
k = 0;
f0 = @(x_prev) (x_prev/2 + 25*x_prev/(1+x_prev.^2));
fk = @(x_prev,k) (x_prev/2 + 25*x_prev/(1+x_prev.^2) + 8*cos(1.2*k));
cov = 10;
propagation = ConditionalNonlinearWithAdditiveGaussian(f0, cov);

y = @(x) (x.^2/20);
cov = 1;
measurement = ConditionalNonlinearWithAdditiveGaussian(y, cov);

% Create filter testing
%sis = BootstrapSISFilter(propagation, @(x,y)propagation.pdf(x,y), @(x,y)measurement.pdf(x,y));
%sis = SIRFilter(propagation, @(x,y)propagation.pdf(x,y), @(x,y)measurement.pdf(x,y));
%uniform_distribution = Uniform(-15, 25);
%sis = ImprovedSIRFilter(propagation, @(x,y)propagation.pdf(x,y), @(x,y)measurement.pdf(x,y), uniform_distribution, 0*0.05);
sis = AdaptiveParticleFilter(propagation, @(x,y)propagation.pdf(x,y), @(x,y)measurement.pdf(x,y));
sis = sis.initKLD(-15, 25, 75, 0.05, 0.05, 10);

% Initialize particles by Importance sampling a uniform distribution center
% along a Gaussian target PDF with mean around the true test state
n_particles = 1000;
x_true = 17;
target_pdf = Gaussian(x_true, 2);
proposal_distribution = Uniform(-15, 25);
[x, weights] = ImportanceSampling_WeightedSamples(@(x)target_pdf.pdf(x), proposal_distribution, n_particles);
sis.particles = x;
sis.weights = weights;

% Plot initial state
figure(1);
subplot(2,1,1);
sis.plot();
vline(x_true, 'r', 'True state');
xlim([-15, 25]);

subplot(2,1,2);
Histogram(-15, 25, 75, sis.particles, sis.weights).plot();
vline(x_true, 'r', 'True state');
xlim([-15, 25]);

%%
k = k + 1;
x_true = fk(x_true, k);
meas = y(x_true);
% These steps below only happens because the propagation model for this
% benchmark changes for each time-step
prop = ConditionalNonlinearWithAdditiveGaussian(@(x)fk(x,k), cov);
sis.propagation_proposal_distribution = prop;
sis.propagation_pdf = @(x,y)prop.pdf(x,y);
%
sis = sis.filter(meas);
subplot(2,1,1);
sis.plot();
vline(x_true, 'r', 'True state');
xlim([-15, 25]);

subplot(2,1,2);
Histogram(-15, 25, 75, sis.particles, sis.weights).plot();
vline(x_true, 'r', 'True state');
xlim([-15, 25]);