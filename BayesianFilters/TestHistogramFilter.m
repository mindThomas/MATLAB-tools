n = 100;
xmin = -5;
xmax = 20;

mu = 0;
sigma = 1;
initial = Gaussian(mu, sigma^2);

state = initial;
propagation = ConditionalGaussian(1.5, 1, 0.5);

figure(1);
subplot(3,2,1);
initial.plot();
xlim([xmin, xmax]); ylim([0 0.4]);
title('Analytically Propagation');
subplot(3,2,3);
propagation.plot();
xlim([xmin, xmax]); ylim([0 1.5]);

% Propagate analytically (since we are using Gaussians)
subplot(3,2,5);
state.plot();
hold on;
for (i = 1:10)
    % Visualize 10 propagations
    state = propagation.join(state).marginalize(2);
    state.plot();    
end
hold off;
xlim([xmin, xmax]); ylim([0 0.4]);

%% Now try again with the Discrete Bayes filter using a Histogram
hf = HistogramFilter(xmin, xmax, n, @(x,y)propagation.pdf(x,y), @(x)0, @(x)initial.pdf(x));
subplot(3,2,2);
hf.plot();
xlim([xmin, xmax]); ylim([0 0.4]*((xmax-xmin)/n));
title('Histogram (Discrete Bayes) Propagation');
subplot(3,2,4);
Histogram(xmin, xmax, n, @(x)propagation.pdf(x,0)).plot();
xlim([xmin, xmax]); ylim([0 1.5]*((xmax-xmin)/n));
subplot(3,2,6);
hf.plot();
hold on;
for (i = 1:10)
    % Visualize 10 propagations
    hf = hf.predict();
    hf.plot();
end
hold off;
xlim([xmin, xmax]); ylim([0 0.4]*((xmax-xmin)/n));