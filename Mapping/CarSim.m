classdef CarSim
    properties
        dt
        % Motion model
        f
        Fx
        Fu
        Fq
        %
        pose % x, y, theta (heading)
        map
        % Lidar variables
        lidar_beam_sensors % vector of azimuth angles for the LiDAR beams
        latest_scan
        
        robot_radius
    end
    
    methods
        function obj = CarSim(dt, map, angle)
            obj.pose = [0; 0; 0];
            lidar_beams = (0:10:359)';
            for (i = 1:length(lidar_beams))
                azimuth = lidar_beams(i);
                obj.lidar_beam_sensors{i} = BeamRangeSensor(deg2rad(azimuth), deg2rad(2), 20, 0.02);
            end
            %obj.lidar_beam_sensors{1} = BeamRangeSensor(deg2rad(angle), deg2rad(2), 20, 0.2);
            
            obj.dt = dt;
            [obj.f, obj.Fx, obj.Fu, obj.Fq] = CoordinatedTurnModel_Discrete_Thrun(dt);
            
            obj.map = map;
            
            obj.robot_radius = 0.5;
        end
        
        function obj = loadMap(obj, img_path)
            img = imread(img_path);
        end
        
        function obj = move(obj, v, omega, varargin)
            % Implement bicycle model for propagation
            % v = translational velocity
            % rho = steering angle
            % omega = angular velocity            
            if (nargin == 4)
                deterministic = varargin{1};
            else
                deterministic = false;
            end

            % Correlation parameters
            alpha1 = 0.01; % velocity contribution into velocity variance 
            alpha2 = 0; % angular velocity contribution into velocity variance 
            alpha3 = 0.05; % velocity contribution into angular velocity variance 
            alpha4 = 10; % angular velocity contribution into angular velocity variance 
            alpha5 = 0; % velocity contribution into angle pertubation variance 
            alpha6 = 10; % angular velocity contribution into angle pertubation variance 
            % Define pertubation variables
            v_var = (alpha1*v^2 + alpha2*omega^2);
            omega_var = (alpha3*v^2 + alpha4*omega^2);
            gamma_var = (alpha5*v^2 + alpha6*omega^2);            
            gamma = 0;
            
            % Add noise
            if (deterministic == false)
                v = v + sqrt(v_var)*randn(1,1);
                omega = omega + sqrt(omega_var)*randn(1,1);
                gamma = sqrt(gamma_var)*randn(1,1);
            end
            
            % Practical fix against zero angular velocity which will cause
            % division by zero in the motion model
            if (abs(omega) < 100*eps)
                omega = 100*eps;
            end
            
            % Construct motion model inputs
            u = [v; omega];
            q = gamma;
            
            obj.pose = obj.f(obj.pose, u, q);
        end
        
        function scan = getLiDAR2D_Deterministic(obj)
            % Do ray-casting on the map to generate 'lidar_beams' number of
            % beams
            scan = zeros(length(obj.lidar_beam_sensors), 1);
            for (i = 1:length(obj.lidar_beam_sensors))
                scan(i) = obj.lidar_beam_sensors{i}.measurement_deterministic(obj.pose, obj.map);
            end                                     
        end
        
        function scan = getLiDAR2D_Stochastic(obj)
            % Call getLiDAR2D_Deterministic() but add stochastic properties
            % according to the range-bearing measurement model from
            % Probabilistic Robotics
        end
        
        function plot(obj)
            % Plot map
            cla;
            imagesc('CData', 1-obj.map.grid(end:-1:1,:), 'XData', [obj.map.pmin(1), obj.map.pmax(1)], 'YData', [obj.map.pmin(2), obj.map.pmax(2)])
            colormap(gray); % Use a gray colormap
            hold on;
    
            % Plot robot pose    
            % Draw circle of robot center
            pos = [obj.pose(1:2)'-obj.robot_radius 2*obj.robot_radius 2*obj.robot_radius]; % [ [x y] width height ]            
            rectangle('Position', pos, 'Curvature', [1 1], 'EdgeColor', 'b', 'LineWidth', 1);
            % Draw heading direction
            d_heading = obj.robot_radius * [cos(obj.pose(3)); sin(obj.pose(3))];
            p(1:2,1) = obj.pose(1:2);
            p(1:2,2) = p(1:2,1) + d_heading;
            plot(p(1,:), p(2,:), 'b-', 'LineWidth', 2);
                        
            % Plot latest Lidar scan
            for (i = 1:length(obj.latest_scan))
                azimuth = obj.pose(3) + obj.lidar_beam_sensors{i}.azimuth_angle_rad;
                d_radius = obj.robot_radius * [cos(azimuth); sin(azimuth)];
                d_ray = obj.latest_scan(i) * [cos(azimuth); sin(azimuth)];
                % Construct LiDAR ray
                p(1:2,1) = obj.pose(1:2) + d_radius;
                p(1:2,2) = obj.pose(1:2) + d_ray;
                line(p(1,:), p(2,:), 'Color', 'r', 'LineWidth', 1);
            end
            
            hold off;
            axis equal; % Make axes grid sizes equal
            xlim([obj.map.pmin(1), obj.map.pmax(1)]);
            ylim([obj.map.pmin(2), obj.map.pmax(2)]);            
        end
        
        function plotParticles(obj, poses)
            plot(poses(1,:), poses(2,:), '.', 'MarkerSize', 10);
            hold on;
            quiver(poses(1,:), poses(2,:), cos(poses(3,:)), sin(poses(3,:)), 1);
            hold off;
            xlim([obj.map.pmin(1), obj.map.pmax(1)]);
            ylim([obj.map.pmin(2), obj.map.pmax(2)]);
            axis equal;
            grid;
        end
    end
end