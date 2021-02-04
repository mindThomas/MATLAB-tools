q0 = [0,0,0];
q1 = [4,4,deg2rad(-10)];

path1 = ReedsShepp.fit(q0, q1, 0.5, true);
path2 = Dubins.fit2(q0, q1, 0.5);

s = 0:0.01:p.length;

samples = zeros(length(s), 3);
for (i = 1:length(s))
    q = path1.get_pose(s(i));
    samples(i,:) = q';
end

path1.params'

figure(1);
clf;
%path2.draw('r');
hold on;
plot(samples(:,1), samples(:,2), 'b');
hold off;
axis equal;