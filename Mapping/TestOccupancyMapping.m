map = OccupancyGrid([0,0], [20,20], [0.1, 0.1]);
map = map.load('example.bmp');

occupancy = OccupancyGrid([0,0], [20,20], [0.1, 0.1]);

car = CarSim(0.1, map, 45);
car.pose = pose;
car.pose(1:2) = [6, 3.5];
car.pose(3) = 0;

%%
for (i = 1:200)
figure(1);
car = car.move(2, 0.5, true);
car = car.captureLiDAR2D_Deterministic();
car.plot();

measurement_likelihood_field = car.lidar_beam_sensors{1}.inverse_measurement_simple(car.latest_scan(1), car.pose, occupancy);
combined_likelihood_field = measurement_likelihood_field;
for (i = 2:length(car.lidar_beam_sensors))
    measurement_likelihood_field = car.lidar_beam_sensors{i}.inverse_measurement_simple(car.latest_scan(i), car.pose, occupancy);
    for (y = 1:size(measurement_likelihood_field,1))
        for (x = 1:size(measurement_likelihood_field,2))
            if (combined_likelihood_field(y,x) == 0.5)
                combined_likelihood_field(y,x) = measurement_likelihood_field(y,x);
            elseif (measurement_likelihood_field(y,x) ~= 0.5)
                combined_likelihood_field(y,x) = max(combined_likelihood_field(y,x), measurement_likelihood_field(y,x));
            end
        end
    end
    
    %figure(2);
    %imshow(1-measurement_likelihood_field);    
end
occupancy = occupancy.RangeBearing_Update(combined_likelihood_field);

figure(3);
%imshow(1 - (map.grid + measurement_likelihood_field));
occupancy.plot();
drawnow;
end

%%
for (i = 1:200)
    car = car.move(0.5, 0.1, false);
    car.latest_scan = car.getLiDAR2D_Deterministic();
    car.plot();
    drawnow;    
end