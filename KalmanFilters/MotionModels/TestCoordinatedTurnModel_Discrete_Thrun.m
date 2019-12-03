%[f, Fx, Fu] = CoordinatedTurnModel_Discrete2(0.1);
[f, Fx, Fu, Fq] = CoordinatedTurnModel_Discrete_Thrun(0.1);

x0 = [0,0,0]';

v_mean = 0.5;
v_var = 0.001;
omega_mean = 0.4;
omega_var = 0.5;

% Correlation parameters
alpha1 = 0.01; % velocity contribution into velocity variance 
alpha2 = 0; % angular velocity contribution into velocity variance 
alpha3 = 0.05; % velocity contribution into angular velocity variance 
alpha4 = 10; % angular velocity contribution into angular velocity variance 
alpha5 = 0; % velocity contribution into angle pertubation variance 
alpha6 = 10; % angular velocity contribution into angle pertubation variance 
% Define pertubation variables
v = Gaussian(v_mean, (alpha1*v_mean^2 + alpha2*omega_mean^2));
omega = Gaussian(omega_mean, (alpha3*v_mean^2 + alpha4*omega_mean^2));
gamma = Gaussian(0, (alpha5*v_mean^2 + alpha6*omega_mean^2));

x = [];
for (i = 1:100)
    u = [v.draw();
         omega.draw()];
    q = gamma.draw();
    x(:,i) = f(x0, u, q);
end

plot(x(1,:), x(2,:), '.', 'MarkerSize', 10);
hold on;
quiver(x(1,:), x(2,:), cos(x(3,:)), sin(x(3,:)), 1);
hold off;
ylim([-0.15, 0.15]);
axis equal;
grid;