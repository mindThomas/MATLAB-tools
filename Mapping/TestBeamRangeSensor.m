grid = OccupancyGrid([0,0], [20,20], [0.1, 0.1]);
grid = grid.load('example.bmp');
grid.plot();
%%
grid.grid(80,:) = 1;
pose = [0; 0; deg2rad(0)];

sensor = BeamRangeSensor(deg2rad(45), deg2rad(3), 20, 0.01);

%invMap = sensor.inverse_measurement(2, pose, grid);
%imshow(1-invMap(end:-1:1,:));

grid.plot();
sensor.ray_tracing_hit(pose, grid)