% Inspired from https://gieseanw.wordpress.com/2012/10/21/a-comprehensive-step-by-step-tutorial-to-computing-dubins-paths/

q0 = [0,0,deg2rad(20)];
q1 = [3, 2, deg2rad(90)];

cmax = 0.5;
R = 1/cmax;

figure(1);
clf;

c1 = q0(1:2)' + R*[-sin(q0(3)); cos(q0(3))];
c3 = q1(1:2)' + R*[-sin(q1(3)); cos(q1(3))];

v = c3 - c1;
D = norm(v);
version = 1;
if (version == 1)
    theta = atan2(v(2), v(1)) - acos(D / (4*R));
elseif (version == 2)
    theta = atan2(v(2), v(1)) + acos(D / (4*R));
end

c2 = c1 + 2*R*[cos(theta); sin(theta)];
u = c2 - c3;
psi = atan2(u(2), u(1))

phi1 = theta + pi/2 - q0(3)
phi2 = theta - psi
phi3 = 3*pi/2 + q1(3) - psi

phi1 = mod(phi1, 2*pi);
phi2 = mod(phi2, 2*pi);
phi3 = mod(phi3, 2*pi);

path = Dubins.from_word(q0, 'LRL', [cmax, phi1, cmax, -phi2, cmax, phi3]);

plot(q0(1)+0.5*[0, cos(q0(3))], q0(2)+0.5*[0, sin(q0(3))]);
hold on;
plot(q1(1)+0.5*[0, cos(q1(3))], q1(2)+0.5*[0, sin(q1(3))]);

plot(c1(1), c1(2), 'bx');
plot(c3(1), c3(2), 'rx');
plot(c2(1), c2(2), 'gx');
circle(c1(1), c1(2), R, 'b');
circle(c3(1), c3(2), R, 'r');
circle(c2(1), c2(2), R, 'g');

path.draw('m')

hold off;

xlim([-5, 5]);
ylim([-5, 5]);
axis equal;

%%
path1 = Dubins.fit(q0, q1, cmax);
path2 = Dubins.fit2(q0, q1, cmax);

path1.params
path2.params