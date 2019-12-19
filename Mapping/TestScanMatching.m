% Load map
map = OccupancyGrid([0,0], [20,20], [0.1, 0.1]);
map = map.load('straight2.bmp');

map_true = OccupancyGrid([0,0], [20,20], [0.1, 0.1]);
map_true = map_true.load('straight3.bmp');
map = map_true;

% Prepare car sim
dt = 0.1;
car = CarSim(dt, map_true, 45);
car.pose = [2, 5.5, 0]';

% First sample
car = car.captureLiDAR2D_Deterministic();
 
% Visualize
figure(1);
car.plot();

figure(2);
plot(car.latest_scan_points(:,1), car.latest_scan_points(:,2), '.');
grid;

%% Move car
v = 0.2;
omega = 0.1;

% Store old scan
prev_scan = car.latest_scan_points;
R = eye(2);
t = [0;0];

%% Move sim car
% Compute odometry transform
[Rdiff, tdiff] = car.odometry_transform(v, omega);
R = R * Rdiff;
t = Rdiff*t + tdiff

% Move sim car
car = car.move(v, omega, false);
% Capture new lidar scan
car = car.captureLiDAR2D_Deterministic();
% Visualize
figure(1);
car.plot();

% Visualize the two scans together
figure(2);
plot(car.latest_scan_points(:,1), car.latest_scan_points(:,2), 'r.');
hold on;
plot(prev_scan(:,1), prev_scan(:,2), 'b.');
hold off;
grid;
axis equal; 

%%
A = car.latest_scan_points';
B = prev_scan';
[R2, t2] = ScanICP(A, B, R, t)

errOdom = car.pose(1:2) - t - [2;5.5]
errICP = car.pose(1:2) - t2 - [2;5.5]