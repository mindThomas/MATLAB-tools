% Load map
map = OccupancyGrid([0,0], [20,20], [0.1, 0.1]);
map = map.load('straight2.bmp');

map_true = OccupancyGrid([0,0], [20,20], [0.1, 0.1]);
map_true = map_true.load('straight3.bmp');
%map = map_true;

% Prepare car sim
dt = 0.1;
car = CarSim(dt, map_true, 45);
car.pose = [2, 5.5, 0]';

% Create Adaptive Particle Filter
localization = AMCL2(map, dt, 30, car.pose);

% Visualization
car.plot();
hold on;
car.plotParticles(localization.particles);
hold off;

%% First sample
car = car.captureLiDAR2D_Deterministic();
%localization = localization.filter(0, 0, car.latest_scan, car.lidar_beam_sensors);
localization = localization.propagate(0, 0);
localization = localization.update(car.latest_scan, car.lidar_beam_sensors);
% 
% car.plot();
% hold on;
% car.plotParticles(localization.particles);
% hold off;

%% Move and sample
v = 0.2;
for (i = 0:100)      
    figure(1);
    for (j = 0:20)
        omega = ((20*i+j) < 500) * 0.2 * cos(2*pi * (20*i+j) / 200);
        %omega = 0;
        
        car = car.move(v, omega, true);
        car = car.captureLiDAR2D_Deterministic();
        localization = localization.propagate(v, omega);

        car.plot();
        hold on;
        car.plotParticles(localization.particles);
        hold off;   
        drawnow;
    end
    
    range = 0:0.01:25;
    p = zeros(length(range),1);
    for (j = 1:length(range))
        p(j) = car.lidar_beam_sensors{5}.inverse_measurement_probability(range(j), car.pose, map_true);
    end
    meas = zeros(length(localization.particles),1);
    meas_p = zeros(length(localization.particles),1);
    for (j = 1:length(localization.particles))
        meas(j) = car.lidar_beam_sensors{5}.measurement_deterministic(localization.particles(:,j), map);
        meas_p(j) = car.lidar_beam_sensors{5}.inverse_measurement_probability(car.latest_scan(5), localization.particles(:,j), map);
    end    
    figure(2);
    plot(range, p);
    hold on;
    stem(meas, meas_p);
    hold off;    
    
    figure(1);
    %car = car.captureLiDAR2D_Deterministic();
    localization = localization.update(car.latest_scan, car.lidar_beam_sensors);

    car.plot();
    hold on;
    car.plotParticles(localization.particles);
    hold off;
    
    range = 0:0.01:25;
    p = zeros(length(range),1);
    for (j = 1:length(range))
        p(j) = car.lidar_beam_sensors{5}.inverse_measurement_probability(range(j), car.pose, map_true);
    end
    meas = zeros(length(localization.particles),1);
    meas_p = zeros(length(localization.particles),1);
    for (j = 1:length(localization.particles))
        meas(j) = car.lidar_beam_sensors{5}.measurement_deterministic(localization.particles(:,j), map);
        meas_p(j) = car.lidar_beam_sensors{5}.inverse_measurement_probability(car.latest_scan(5), localization.particles(:,j), map);
    end    
    figure(2);
    plot(range, p);
    hold on;
    stem(meas, meas_p);
    hold off;       
    
    pause(0.2);    
end