f = @(x, theta) theta(1) + theta(2)*x - theta(3)*cos(theta(4) * x);
Jf = @(x, theta) [ones(length(x), 1),...
                  x,...                 
                  -cos(theta(4) * x),...
                  theta(3)*x.*sin(theta(4) * x)];    

% f = @(x, theta) theta(1) + theta(2)*x - cos(theta(3) * x);
% Jf = @(x, theta) [ones(length(x), 1),...
%                   x,...                                   
%                   x.*sin(theta(3) * x)];                  
              
theta_true = [1, -0.5, 10, 0.6]';
%theta_true = [1, -0.5, 1]';

x = (0:0.2:50)';
y_true = f(x, theta_true);
sigma2 = 2;
y_meas = y_true + sqrt(sigma2)*randn(length(y_true),1);

figure(1);
fplot(@(x)f(x,theta_true), [0, 50]);
hold on;
plot(x, y_meas, '*');
hold off;

%%
% Initial parameter guess by linear fitting
%a = mean(diff(y_meas) ./ diff(x));
%b = mean(y_meas);
A = [x.*x, x, ones(length(x), 1)];
b = y_meas;
sol = A \ b;
a = sol(2);
b = sol(3);

%theta0 = [b, a, 0, 0];
theta0 = [sol(3), sol(2), 0.1, 0.6];
theta0 = [0.2, -0.5, 0.1, 4.5*2*pi/(max(x)-min(x))]';
%theta0 = [0.3, -0.5, 0.8]';
y_guess = f(x, theta0);

figure(2);
fplot(@(x)f(x,theta_true), [0, 50]);
hold on;
plot(x, y_guess, '*');
hold off;

theta_sol = theta0;

% Fit
%theta_sol = NonlinearOptimizers.gradient_descent_data_fitting(theta_sol, @(theta)f(x,theta), @(theta)Jf(x,theta), y_meas, 0.01, 0.00001)
%theta_sol = NonlinearOptimizers.gauss_newton_data_fitting(theta_sol, @(theta)f(x,theta), @(theta)Jf(x,theta), y_meas)
theta_sol = NonlinearOptimizers.levenberg_marquardt(theta_sol, @(theta)f(x,theta), @(theta)Jf(x,theta), y_meas, 1e-8)
y_sol = f(x, theta_sol);

figure(2);
fplot(@(x)f(x,theta_true), [0, 50]);
hold on;
plot(x, y_sol, '*');
hold off;