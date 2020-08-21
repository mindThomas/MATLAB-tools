%% Test basis function
spline = Bspline_simple_closed(0:9, 3);
%spline.knots = 0:9;
%spline.n = length(spline.knots);

%u = min(spline.knots):0.02:max(spline.knots);
u = 0:0.02:10;
B0 = zeros(size(u));
B1 = zeros(size(u));
for (i = 1:length(u))
    %B0(i) = spline.basis1(  7, u(i));
    B1(i) = spline.basis(3, 5, u(i));       
end

figure(1);
plot(u, B0);
hold on;
plot(u, B1);
hold off;

%% Comparison of basis function to Gaussian
degree = 3;
mu = 3;
sigma = 0.6;

spline = Bspline(0:9, degree);
%spline.knots = 0:9;
%spline.n = length(spline.knots);


%u = min(spline.knots):0.02:max(spline.knots);
u = 0:0.02:10;
B1 = zeros(size(u));
B2 = zeros(size(u));
for (i = 1:length(u))    
    B1(i) = spline.basis(degree, 4, u(i));
    B2(i) = 1/(sigma*sqrt(2*pi)) * exp(-1/2 * ((u(i) - mu)/sigma)^2);
end

figure(1);
plot(u, B1);
hold on;
plot(u, B2);
hold off;
xlim([0, 6]);
legend('3-order Bspline basis', 'Gaussian, \sigma=0.6');
vline([mu-3*sigma, mu+3*sigma], {'r--', 'r--'}, {'-3*\sigma', '3*\sigma'});
title('Basis function comparison to a Gaussian distribution');

%% Test 1D curve
control_points = [0.0, 0.0, 2.5, 0.5, 10.0, 10.0, 10, 0.0 ,0.0];
spline = Bspline(control_points, 3);

figure(1);
plot(spline.knots, spline.control_points, 'o');

u = 0:0.01:(length(control_points)-1);
y = zeros(size(u));
for (i = 1:length(u))
   y(i) = spline.eval(u(i));
end

hold on;
plot(u, y);
hold off;
xlim([0, (length(control_points)-1)]);

%% Test 2D curve
type = 1;
s = 1;

if type == 1
control_points = [0, 0, 0, 0, 4, 4, 4, 4;
                  0, 0, 4, 4, 4, 4, 0, 0];
elseif type == 2              
control_points = [0, 0, 4, 4;
                  0, 4, 4, 0];              
elseif type == 3              
control_points = [0, 0, (2-s*2), 2, s*4, 4, s*4, 2;
                  0, 2, (2+s*2), 4, s*4, 2, s*0, 0];                            
end

weights = ones(1,length(control_points));
weights(1) = 6;
weights(2) = 1;
weights(3) = 6;
weights(4) = 1;
weights(5) = 6;
weights(6) = 1;
weights(7) = 6;
weights(8) = 1;

spline = Bspline_simple_closed(control_points, 3);
%spline = NURB_Curve(control_points, weights, 3);

figure(1);
plot(spline.control_points(1,:), spline.control_points(2,:), 'o');
axis equal;

u = 0:0.01:(length(control_points));
y = zeros(2, length(u));
dy = zeros(2, length(u));
heading = zeros(1, length(u));
dheading = zeros(1, length(u));
for (i = 1:length(u))
   y(:,i) = spline.eval(u(i));
   %dy(:,i) = spline.deval(u(i));   
   %heading(:,i) = spline.heading(u(i));   
   %dheading(:,i) = spline.dheading(u(i));   
end
heading = unwrap(heading);

hold on;
plot(y(1,:), y(2,:));
quiver(y(1,1:10:end), y(2,1:10:end), dy(1,1:10:end), dy(2,1:10:end), 0.2);
hold off;
xlim([-0.5, 4.5]);
ylim([-0.5, 4.5]);

%%
dp_du = diff(y') ./ diff(u');
heading_numerical = atan2(dp_du(:,2), dp_du(:,1));
heading_numerical = unwrap(heading_numerical);

dheading_numerical = diff(heading) ./ diff(u);

figure(2);
subplot(2,1,1);
plot(u, rad2deg(heading));
hold on;
plot((u(1:end-1)+u(2:end))/2, rad2deg(heading_numerical));
hold off;
subplot(2,1,2);
plot(u, rad2deg(dheading));
hold on;
plot((u(1:end-1)+u(2:end))/2, rad2deg(dheading_numerical));
hold off;