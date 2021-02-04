clf;

t = ((spline.n-spline.original_n)/2):0.02:(spline.original_n-1+(spline.n-spline.original_n)/2);
evaluation_points = zeros(2, length(t));
for (i = 1:length(t))    
    evaluation_points(:,i) = spline.evaluate(t(i));
end

figure(1);
plot(spline.control_points(1,:), spline.control_points(2,:), 'bo');
hold on;
plot(evaluation_points(1,:), evaluation_points(2,:), 'g');
hold off;
axis equal;
drawnow;

p = [7.5; -1];
t = 12;

%%
q = spline.evaluate(t);
v = p - q

hold on;
plot(p(1), p(2), 'rx');
plot(q(1), q(2), 'bx');
hold off;

kappa = spline.curvature(t)
R = 1/kappa;

df_dt = spline.deval_dt(spline.p, spline.control_points, t);
df_dt = df_dt / norm(df_dt);

df_dt_lat = [0,1;-1,0] * df_dt;

lon_err = df_dt' * v
lat_err = ([0,1;-1,0] * df_dt)' * v

hold on;
plot(q(1)+[0, df_dt(1)], q(2)+[0, df_dt(2)]);
plot(q(1)+[0, df_dt_lat(1)], q(2)+[0, df_dt_lat(2)]);
hold off;

if (abs(lat_err) < abs(R))
    theta = atan(lon_err / (R+lat_err));
else
    theta = sign(lon_err) * deg2rad(90);
end
arc_curve_length_opt = R * theta;
stretch = spline.cartesian_stretch(spline.p, spline.control_points, t)
t_opt = t + arc_curve_length_opt / stretch

q_opt = spline.evaluate(t_opt);
hold on;
plot(q_opt(1), q_opt(2), 'mx');
hold off;

%
t = t_opt;