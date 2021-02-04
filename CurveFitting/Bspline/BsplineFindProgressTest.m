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

while (1)
    [x,y] = ginput(1);
    p = [x;y];
    
    t = spline.find_t(p)
    t2 = Bspline_uniform.find_locations5(spline.p, p(1), p(2), t, spline.control_points(1,:), spline.control_points(2,:))
    
    q = spline.evaluate(t);
    q2 = spline.evaluate(t2);

    hold on;
    plot(p(1), p(2), 'rx');
    plot(q(1), q(2), 'bx');
    plot(q2(1), q2(2), 'mo');
    hold off;
end