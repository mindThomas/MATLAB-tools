Mu = [2;2];
Cov = [2, -1.8
      -1.8, 2];

x = Gaussian(Mu, Cov);
f = @(x) [sqrt(x(1)^2+x(2)^2); atan2(x(2),x(1))];

nRealizations = 100000;
realizations = zeros(size(Mu,1), nRealizations);
propagated = zeros(size(Mu,1), nRealizations);
for (i = 1:nRealizations)
    realizations(:,i) = x.draw();
    propagated(:,i) = f(realizations(:,i));
end

figure(1);
subplot(2,1,1);
plot(realizations(1,:), realizations(2,:), '*');
hold on;
x.plotSigmaContour(3);
hold off;
xlim([-4 8]);
ylim([-4 8]);
axis equal;
title('Input');

subplot(2,1,2);
plot(propagated(1,:), propagated(2,:), '*');
xlim([0 10]);
ylim([-1.5 3]);
hold on;

[Mu2, Cov2] = MomentPropagation.Propagate_Linearized(Mu, Cov, f)
y = Gaussian(Mu2, Cov2);
y.plotSigmaContour(3, 'k-');

[Mu2, Cov2] = MomentPropagation.Propagate_Unscented(Mu, Cov, f)
y = Gaussian(Mu2, Cov2);
y.plotSigmaContour(3, 'r-');

[Mu2, Cov2] = MomentPropagation.Propagate_Cubature(Mu, Cov, f)
y = Gaussian(Mu2, Cov2);
y.plotSigmaContour(3, 'g-');

Mu2 = mean(propagated, 2);
Cov2 = cov(propagated')';
y = Gaussian(Mu2, Cov2);
y.plotSigmaContour(3, 'c-');
hold off;

title('Output Comparisons');
legend('Realizations', 'Linearized', 'Unscented', 'Cubature', 'Sample based');

%%
figure(2);
[Mu2, Cov2, sigmaPointsIn, sigmaPointsOut] = MomentPropagation.Propagate_Unscented(Mu, Cov, f)
y = Gaussian(Mu2, Cov2);

subplot(2,2,1);
plot(realizations(1,:), realizations(2,:), '*');
hold on;
x.plotSigmaContour(3);
plot(sigmaPointsIn(1,1), sigmaPointsIn(2,1), 'go');
plot(sigmaPointsIn(1,2:end), sigmaPointsIn(2,2:end), 'ro');
hold off;
xlim([-4 8]);
ylim([-4 8]);
axis equal;
title('Unscented Transform');

subplot(2,2,3);
plot(propagated(1,:), propagated(2,:), '*');
hold on;
y.plotSigmaContour(3, 'r-');
plot(sigmaPointsOut(1,1), sigmaPointsOut(2,1), 'go');
plot(sigmaPointsOut(1,2:end), sigmaPointsOut(2,2:end), 'ro');
hold off;
xlim([0 10]);
ylim([-1.5 3]);

[Mu2, Cov2, sigmaPointsIn, sigmaPointsOut] = MomentPropagation.Propagate_Cubature(Mu, Cov, f)
y = Gaussian(Mu2, Cov2);

subplot(2,2,2);
plot(realizations(1,:), realizations(2,:), '*');
hold on;
x.plotSigmaContour(3);
plot(sigmaPointsIn(1,:), sigmaPointsIn(2,:), 'ro');
hold off;
xlim([-4 8]);
ylim([-4 8]);
axis equal;
title('Cubature rule');

subplot(2,2,4);
plot(propagated(1,:), propagated(2,:), '*');
hold on;
y.plotSigmaContour(3, 'r-');
plot(sigmaPointsOut(1,:), sigmaPointsOut(2,:), 'ro');
hold off;
xlim([0 10]);
ylim([-1.5 3]);