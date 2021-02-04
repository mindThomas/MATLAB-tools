figure(3);
clf;

q0 = [0;0;0];
%p1 = [ 0.5,deg2rad(100), 2, 0.2,deg2rad(-45) ];
p1 = [ 0.1,deg2rad(0), 2, 0.2,deg2rad(90) ];
Dubins.draw_csc(q0,  p1(1),p1(2),p1(3),p1(4),p1(5), 'b');
[x,y,theta] = Dubins.dubins_csc(q0, p1(1),p1(2),p1(3),p1(4),p1(5))
[x,y,theta] = Dubins.dubins_csc_compact(q0, p1(1),p1(2),p1(3),p1(4),p1(5))

dubins_length1 = Dubins.arc_length(p1(1),p1(2)) + p1(3) + Dubins.arc_length(p1(4),p1(5));

q1 = [x;y;theta];

%[p2, dubins_length2] = Dubins.find_parameters_CSC(q0, q1, 0.5);
%Dubins.draw_csc(q0,  p2(1),p2(2),p2(3),p2(4),p2(5), 'r');

[p2, dubins_length2] = Dubins.find_parameters_CCC(q0, q1, 0.5, false);
Dubins.draw_ccc(q0,  p2(1),p2(2),p2(3),p2(4),p2(5),p2(6), 'r');


%%
b = Dubins.fit(q0, q1);
b.draw('g');

poses = b.get_pose(0:0.5:b.length);
hold on;
plot(poses(1,:), poses(2,:), 'xb');
hold off;

axis equal;

dubins_length1
dubins_length2