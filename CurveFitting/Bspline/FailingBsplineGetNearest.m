points = [[
                -15.140398,
                4.850069,
                1.008874
            ],...
        [
                -9.251683,
                6.831848,
                0.608689
            ],...
        [
                4.570328,
                5.847807,
                -0.173246
            ],...
        [
                -0.984834,
                12.777185,
                1.06663
            ],...
        [
                0.365669,
                20.619912,
                0.670122
            ],...
        [
                -8.591163,
                22.023434,
                3.141592653589793
            ],...
        [
                -19.141724,
                22.246333,
                -2.303612
            ],...
        [
                -22.188013,
                11.547172,
                -2.096863
            ],...
        [
                -12.974847,
                10.544126,
                0.0
            ] ];    

points = points(1:2,:);
control_points = [points(:,1), points(:,1), points, points(:,end), points(:,end)];
        
bspline = Bspline_uniform();
bspline.p = 2;
bspline.control_points = control_points;
bspline.n = size(control_points, 2);
bspline.original_control_points = points;
bspline.original_n = size(points, 2);

figure(1);
bspline.plot();

%goal = [3.621404; 16.87674];
%goal = [1.60241; 7.69819];
%goal = [1.258658; 11.859199];
%goal = [-18.763354; 15.367815];
%goal = [-15.128186; 15.710371];
goal = [-16.983156; 15.848195];

hold on;
plot(goal(1), goal(2), 'ro');
hold off;

t_goal = bspline.find_t(goal)
goal_evaluated = bspline.evaluate(t_goal);
hold on;
plot(goal_evaluated(1), goal_evaluated(2), 'rx');
hold off;