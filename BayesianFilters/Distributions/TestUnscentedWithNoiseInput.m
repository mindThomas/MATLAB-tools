Mu = [2;2];
Cov = [2, -1.8
      -1.8, 2];

x = Gaussian(Mu, Cov);
r = Gaussian(0, 0.5);
f = @(x,r) [sqrt(x(1)^2+x(2)^2)+r(1); atan2(x(2),x(1))];

nRealizations = 100000;
realizations = zeros(size(Mu,1), nRealizations);
propagated = zeros(size(Mu,1), nRealizations);
for (i = 1:nRealizations)
    realizations(:,i) = x.draw();
    propagated(:,i) = f(realizations(:,i), r.draw());
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


%%
r0 = 0;
Mu2 = f(Mu,r0);
[dfdx, dfdr] = numjacobian2_real(f, Mu, r0);
Cov2 = dfdx * x.Cov * dfdx' + dfdr * r.Cov * dfdr';
y = Gaussian(Mu2, Cov2);
y.plotSigmaContour(3, 'y-');

Mu_prime = [Mu; r.Mu];
Cov_prime = [x.Cov, zeros(size(x.Cov,1), size(r.Cov,2));
             zeros(size(r.Cov,1), size(x.Cov,2)), r.Cov];
nx = size(x.Mu,1);
[Mu2, Cov2] = MomentPropagation.Propagate_Linearized(Mu_prime, Cov_prime, @(x)f(x(1:nx), x(nx+1:end)))
y = Gaussian(Mu2, Cov2);
y.plotSigmaContour(3, 'k-');

[Mu2, Cov2] = MomentPropagation.Propagate_Unscented(Mu_prime, Cov_prime, @(x)f(x(1:nx), x(nx+1:end)))
y = Gaussian(Mu2, Cov2);
y.plotSigmaContour(3, 'r-');

[Mu2, Cov2] = MomentPropagation.Propagate_Cubature(Mu_prime, Cov_prime, @(x)f(x(1:nx), x(nx+1:end)))
y = Gaussian(Mu2, Cov2);
y.plotSigmaContour(3, 'g-');

Mu2 = mean(propagated, 2);
Cov2 = cov(propagated')';
y = Gaussian(Mu2, Cov2);
y.plotSigmaContour(3, 'c-');
hold off;

title('Output Comparisons');
legend('Realizations', 'Linearized', 'Unscented', 'Cubature', 'Sample based');
