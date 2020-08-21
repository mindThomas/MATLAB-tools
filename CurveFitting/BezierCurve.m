%% See https://en.wikipedia.org/wiki/B%C3%A9zier_curve

%% 1D Quadratic Bezier curve example
% Define the 3 control points
P = [2, 0, 1]

B = @(t) ((1-t).^2 * P(1) + 2*(1-t).*t * P(2) + t.^2 * P(3));

figure(1);
plot([0, 0.5, 1], P, 'o');
hold on;
line([0, 0.5], P(1:2), 'LineStyle', '--');
line([0.5, 1], P(2:3), 'LineStyle', '--');
fplot(B, [0,1]);
hold off;

%% 2D Quadratic Bezier curve example
% Define the 3 control points
P = [2, 0.5;
     0, 0;
     0.5, 2.5]

B = @(t) ((1-t).^2 .* P(1,:)' + 2*(1-t).*t .* P(2,:)' + t.^2 .* P(3,:)');
% Polynomial formulation
%B = @(t) (P(1,:)' - 2*(P(1,:)' - P(2,:)') .* t + ((P(3,:)' - P(2,:)') + (P(1,:)' - P(2,:)')) .* t.^2);

figure(1);
plot(P(:,1), P(:,2), 'o');
hold on;
line(P(1:2, 1), P(1:2, 2), 'LineStyle', '--');
line(P(2:3, 1), P(2:3, 2), 'LineStyle', '--');
t = 0:0.01:1;

t = [0; 1];
for (n = 1:100)
    t = [t;rand];
end
t = sort(t)';

curve = B(t);
plot(curve(1,:), curve(2,:));
hold off;

grid;
axis equal;
xlim([-1, 3]);
ylim([-1, 3]);


%% 2D Cubic Bezier curve example
% Define the 4 control points
P = [2, 0.5;
     0, 0;
     -1, 1;
     0.5, 2.5]

B = @(t) ((1-t).^3 .* P(1,:)' + 3*(1-t).^2.*t .* P(2,:)' + 3*(1-t).*t.^2 .* P(3,:)' + t.^3 .* P(4,:)');

figure(1);
plot(P(:,1), P(:,2), 'o');
hold on;
line(P(1:2, 1), P(1:2, 2), 'LineStyle', '--');
line(P(3:4, 1), P(3:4, 2), 'LineStyle', '--');
t = 0:0.01:1;
curve = B(t);
plot(curve(1,:), curve(2,:));
hold off;

grid;
axis equal;
xlim([-1, 3]);
ylim([-1, 3]);


%% Fit 2D Quadratic Bezier curve example
% Define the 3 control points
P = [2, 0.5;
     0, 0;
     0.5, 2.5]

B = @(P0,P1,P2,t) ((1-t).^2 .* P0 + 2*(1-t).*t .* P1 + t.^2 .* P2);

% Generate samples from random locations from the curve
t = [0; 1];
for (n = 1:100)
    t = [t;rand];
end
t = sort(t)';
curve = B(P(1,:)', P(2,:)', P(3,:)', t);

% Add random noise to the points
x = curve(1,:)' + 0.01*randn(length(t),1);
y = curve(2,:)' + 0.01*randn(length(t),1);

% Fit Bezier curve
P = fitBezierCurve(x, y);

% Evaluate fitted curve at linearly spaced intervals
t = linspace(0, 1, length(t));
curve_fitted = B(P(1,:)', P(2,:)', P(3,:)', t);

% Plot results
figure(1);    
plot(x, y, '.'); % data points
hold on;
plot(curve(1,:), curve(2,:), '--'); % original curve    

line(P(1:2, 1), P(1:2, 2), 'LineStyle', '--');
line(P(2:3, 1), P(2:3, 2), 'LineStyle', '--');
plot(P(:,1), P(:,2), 'o');

plot(curve_fitted(1,:), curve_fitted(2,:));
hold off;

grid;
axis equal;
xlim([-1, 3]);
ylim([-1, 3]);