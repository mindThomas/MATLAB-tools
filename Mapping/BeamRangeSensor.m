classdef BeamRangeSensor
    % Implements the Beam range sensor fom "Probabilistic Robotics"
    % have four different terms as part of the measurement distribution:
    % p_hit, p_short, p_max, p_rand
    properties
       azimuth_angle_rad
       beam_width_rad      
       z_max
       var_hit % range measurement variance when hitting the object
       lambda_short % exponential rate for measurement likelihood to increase towards the true value
       % Weighting factors for the different distributions
       w_hit
       w_short
       w_max
       w_rand
    end
    
    methods

        function obj = BeamRangeSensor(azimuth, beam_width, max_range, var_hit)            
            obj.azimuth_angle_rad = azimuth;
            obj.beam_width_rad = beam_width;            
            obj.z_max = max_range;
            obj.var_hit = var_hit; 
            
            obj.lambda_short = 0.5;
            obj.w_hit = 1.0;
            obj.w_short = 0.75;
            obj.w_max = 0.5;
            obj.w_rand = 0*1.0;
            
            % Normalize weights
            eta = obj.w_hit + obj.w_short + obj.w_max + obj.w_rand;
            obj.w_hit = obj.w_hit / eta;
            obj.w_short = obj.w_short / eta;
            obj.w_max = obj.w_max / eta;
            obj.w_rand = obj.w_rand / eta;
        end
        
        function z = measurement(obj, pose, map)
            % Call measurement_deterministic() but add stochastic properties
            % according to the range-bearing measurement model from
            % Probabilistic Robotics
            z_deterministic = obj.measurement_deterministic(pose, map);
        end
                 
        function z = measurement_deterministic(obj, pose, map)
            z = obj.ray_tracing_hit(pose, map);
        end               
        
        function [cell_value, cell_found] = inverse_measurement_hit(obj, range_measurement, pose, map)
            % Do a simple lookup in the map                        
            % Compute the location of the captured cell            
            x_z = pose(1) + range_measurement * cos(pose(3) + obj.azimuth_angle_rad);
            y_z = pose(2) + range_measurement * sin(pose(3) + obj.azimuth_angle_rad);
            
            if (x_z >= map.pmin(1) && x_z <= map.pmax(1) && ...
                y_z >= map.pmin(2) && y_z <= map.pmax(2))
                [row,col] = map.getRowCol(x_z, y_z);
                cell_found = true;   
                cell_value = map.grid(row,col);
                return;
            end
            
            % Cell is outside the map
            cell_found = false;
            cell_value = 0.5;                                      
        end 
        
        function p = inverse_measurement_probability(obj, range_measurement, pose, map)
            range = obj.ray_tracing_hit(pose, map);
            
            p = (...
                obj.w_hit   * obj.p_hit(range_measurement, range) ...
              + obj.w_short * obj.p_short(range_measurement, range) ...
              + obj.w_max   * obj.p_max(range_measurement, range) ...
              + obj.w_rand  * obj.p_rand(range_measurement, range) ...
              );            
        end                 
        
        function probability_map = inverse_measurement(obj, range_measurement, pose, occupancy_map)
            % Generate an inverse measurement model probability map
            % https://www.coursera.org/lecture/motion-planning-self-driving-cars/lesson-2-populating-occupancy-grids-from-lidar-scan-data-part-2-VcH67
            % Occupied = 1.0
            % No information/unknown = 0.5
            % Free = 0.0
            mask = zeros(size(occupancy_map.grid));
            probability_map = zeros(size(occupancy_map.grid));           
            
            % Loop over all cells and compute probability
            for (i = 1:size(probability_map, 1))
                for (j = 1:size(probability_map, 2))
                    [x_m, y_m] = occupancy_map.getXY(i,j);
                    probability_map(i,j) = obj.inverse_measurement_xy(range_measurement, pose, x_m, y_m);
                    if (probability_map(i,j) > 0)
                        mask(i,j) = 1;
                    end
                end
            end
            
            % Normalize map
            %probability_map = probability_map / max(max(probability_map));
            
            % Add outside of FoV unknown probability of 0.5
            probability_map = probability_map + 0.5*(1-mask);
        end             
        
        function probability_map = inverse_measurement_simple(obj, range_measurement, pose, occupancy_map)
            % Generate an inverse measurement model probability map
            % See Table 9.2 in Probabilistic Robotics
            %
            % Occupied = 1.0
            % No information/unknown = 0.5
            % Free = 0.0            
            probability_map = zeros(size(occupancy_map.grid));           
            alpha = sqrt(obj.var_hit);
            
            % Loop over all cells and compute probability
            for (i = 1:size(probability_map, 1))
                for (j = 1:size(probability_map, 2))
                    [x_m, y_m] = occupancy_map.getXY(i,j);
                    range = sqrt( (x_m - pose(1))^2 + (y_m - pose(2))^2 );
                    azimuth = atan2(y_m - pose(2), x_m - pose(1)) - pose(3);
                    if (azimuth < 0)
                        azimuth = azimuth + 2*pi;
                    end                    
                    
                    if ((range > min(obj.z_max, range_measurement+alpha/2)) || (abs(azimuth - obj.azimuth_angle_rad) > obj.beam_width_rad/2))
                        probability_map(i,j) = 0.5; % unknown (outside FOV)
                    elseif ((range_measurement < obj.z_max) && (abs(range - range_measurement) < alpha/2))
                        probability_map(i,j) = 1.0; % occupied
                    elseif (range <= range_measurement)
                        probability_map(i,j) = 0.0; % free
                    else
                        probability_map(i,j) = 0.5;
                    end                    
                end
            end
        end           
        
        function p = inverse_measurement_xy(obj, range_measurement, pose, x_m, y_m)
            % See also Table 9.2 in Probabilistic Robotics
            % Compute the probability of a measurement given a map cell
            % location, (x,y)
            range = sqrt( (x_m - pose(1))^2 + (y_m - pose(2))^2 );
            % consider to take the size of the grid cell into account and
            % the angle onto which the grid cell is observed, such that
            % sqrt(2)/2 of the resolution is substracted from the range if
            % the cell is observed at an angle of 45 degrees since the
            % detection will hit the corner of the cell
            azimuth = atan2(y_m - pose(2), x_m - pose(1)) - pose(3);
            if (azimuth < 0)
                azimuth = azimuth + 2*pi;
            end
            
            if (azimuth >= obj.azimuth_angle_rad-obj.beam_width_rad/2 && azimuth <= obj.azimuth_angle_rad+obj.beam_width_rad/2)
                %eta_beamwidth = 1/obj.beam_width_rad;
                p = (...
                    obj.w_hit   * obj.p_hit(range_measurement, range) ...
                  + obj.w_short * obj.p_short(range_measurement, range) ...
                  + obj.w_max   * obj.p_max(range_measurement, range) ...
                  + obj.w_rand  * obj.p_rand(range_measurement, range) ...
                  );
            else
                p = 0;
            end
        end
        
        function p = measurement_likelihood_field(obj, range_measurement, pose, map)
            % Compute the likelihood of a measurement using the likelihood
            % field approach of finding the nearest occupied cell to the measurement
            x_z = pose(1) + range_measurement * cos(pose(3) + obj.azimuth_angle_rad);
            y_z = pose(2) + range_measurement * sin(pose(3) + obj.azimuth_angle_rad);
            idx = find(map.grid > 0); % find occupied cells
            distances = zeros(size(idx));
            
            % Loop over all occupied cells and compute distance
            for (i = 1:length(idx))
                [x_m, y_m] = map.getXYfromIndex(idx(i));
                distances(i) = sqrt( (x_z - x_m)^2 + (y_z - y_m)^2 );
            end
            
            % Find the closest occupied cell
            [minDist, minDistIdx] = min(distances);
            
            [x_m, y_m] = map.getXYfromIndex(minDistIdx);
            range = sqrt( (x_m - pose(1))^2 + (y_m - pose(2))^2 );
            
            % Compute probability of measurement given the detection of
            % this cell plus the probability of a random detection
            p = obj.w_hit * obj.p_hit(range_measurement, range) ...
              + obj.w_rand * obj.p_rand(range_measurement, range);
        end
                
        function z = ray_tracing_hit(obj, pose, map)
            % Bresenham's line tracing algorithm
            % https://www.coursera.org/lecture/motion-planning-self-driving-cars/lesson-2-populating-occupancy-grids-from-lidar-scan-data-part-2-VcH67
            % From https://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm
            % For other Ray tracing algorithms see also:
            % RangeLibc @ https://github.com/kctess5/range_libc
            [y0, x0] = map.getRowCol(pose(1), pose(2));
            x1 = pose(1) + obj.z_max * cos(pose(3) + obj.azimuth_angle_rad);
            y1 = pose(2) + obj.z_max * sin(pose(3) + obj.azimuth_angle_rad);
            [y1, x1] = map.getRowCol(x1, y1);            
            
            deltax = x1 - x0;
            deltay = y1 - y0;           
                        
            if (abs(deltax) > abs(deltay))
                deltaerr = abs(deltay / deltax); % Assume deltax != 0 (line is not vertical),
                % note that this division needs to be done in a way that preserves the fractional part
                error = 0.0; % No error at start
                y = y0;
                dir = -1 + 2*(x1 > x0);
                for (x = x0:dir:x1)
                    if (x < 1 || x > size(map.grid,2) || y < 1 || y > size(map.grid,1))
                        break;
                    end 
                    hit = map.grid(y,x);                   
                    
                    if (hit > 0.5)
                        % Convert map location back into position                        
                        [x_cell,y_cell] = map.getXY(y,x);
                                                                        
                        min_x = x_cell - pose(1) - map.resolution(1)/2;
                        max_x = x_cell - pose(1) + map.resolution(1)/2;
                        min_y = y_cell - pose(2) - map.resolution(2)/2;
                        max_y = y_cell - pose(2) + map.resolution(2)/2;                        
                        
                        if (abs(min_x) < abs(max_x))
                            closest_x = min_x;
                        else
                            closest_x = max_x;
                        end
                        if (abs(min_y) < abs(max_y))
                            closest_y = min_y;
                        else
                            closest_y = max_y;
                        end                        
                        
                        % Find range where the cell is hit
                        angle = pose(3) + obj.azimuth_angle_rad;                        
                        
                        % OBS! Handle situation where
                        % angle = 0, 90, 180, 270
                        if (angle == 0 || angle == pi)
                            z = abs(closest_x);
                            return;
                        elseif (angle == pi/2 || angle == 3*pi/2)
                            z = abs(closest_y);
                            return;
                        end
                        
                        ratio = sin(angle) / cos(angle);
                        % y = ratio * x
                        % x = y / ratio
                        
                        edge1_y = closest_y;
                        edge1_x = edge1_y / ratio;
                        range1 = sqrt(edge1_x^2 + edge1_y^2);                        
                        if (abs(edge1_x) >= abs(closest_x) && edge1_x >= min_x && edge1_x <= max_x)
                            z = range1;
                            return
                        end
                        
                        edge2_x = closest_x;
                        edge2_y = ratio * edge2_x; 
                        range2 = sqrt(edge2_x^2 + edge2_y^2);
                        if (abs(edge2_y) >= abs(closest_y) && edge2_y >= min_y && edge2_y <= max_y)
                            z = range2;
                            return
                        end                                                                                                
                        
                        % Something is wrong, return default range
                        z = sqrt( (x_cell - pose(1))^2 + (y_cell - pose(2))^2 );
                        return;
                    end
                    error = error + deltaerr;
                    if (error >= 0.5)
                        y = y + sign(deltay) * 1;
                        error = error - 1.0;
                    end
                end
            else
                deltaerr = abs(deltax / deltay); % Assume deltax != 0 (line is not vertical),
                % note that this division needs to be done in a way that preserves the fractional part
                error = 0.0; % No error at start
                x = x0;
                dir = -1 + 2*(y1 > y0);
                for (y = y0:dir:y1)
                    if (x < 1 || x > size(map.grid,2) || y < 1 || y > size(map.grid,1))
                        break;
                    end
                    hit = map.grid(y,x);                   
                    
                    if (hit > 0.5)
                        % Convert map location back into position                        
                        [x_cell,y_cell] = map.getXY(y,x);
                                                                        
                        min_x = x_cell - pose(1) - map.resolution(1)/2;
                        max_x = x_cell - pose(1) + map.resolution(1)/2;
                        min_y = y_cell - pose(2) - map.resolution(2)/2;
                        max_y = y_cell - pose(2) + map.resolution(2)/2;                        
                        
                        if (abs(min_x) < abs(max_x))
                            closest_x = min_x;
                        else
                            closest_x = max_x;
                        end
                        if (abs(min_y) < abs(max_y))
                            closest_y = min_y;
                        else
                            closest_y = max_y;
                        end                        
                        
                        % Find range where the cell is hit
                        angle = pose(3) + obj.azimuth_angle_rad;                        
                        
                        % OBS! Handle situation where
                        % angle = 0, 90, 180, 270
                        if (angle == 0 || angle == pi)
                            z = abs(closest_x);
                            return;
                        elseif (angle == pi/2 || angle == 3*pi/2)
                            z = abs(closest_y);
                            return;
                        end
                        
                        ratio = sin(angle) / cos(angle);
                        % y = ratio * x
                        % x = y / ratio
                        
                        edge1_y = closest_y;
                        edge1_x = edge1_y / ratio;
                        range1 = sqrt(edge1_x^2 + edge1_y^2);                        
                        if (abs(edge1_x) >= abs(closest_x) && edge1_x >= min_x && edge1_x <= max_x)
                            z = range1;
                            return
                        end
                        
                        edge2_x = closest_x;
                        edge2_y = ratio * edge2_x; 
                        range2 = sqrt(edge2_x^2 + edge2_y^2);
                        if (abs(edge2_y) >= abs(closest_y) && edge2_y >= min_y && edge2_y <= max_y)
                            z = range2;
                            return
                        end                                                                                                
                        
                        % Something is wrong, return default range
                        z = sqrt( (x_cell - pose(1))^2 + (y_cell - pose(2))^2 );
                        return;
                    end
                    error = error + deltaerr;
                    if (error >= 0.5)
                        x = x + sign(deltax) * 1;
                        error = error - 1.0;
                    end
                end
            end
            
            z = obj.z_max;
        end
    end
    
    methods (Access = private)
        % Correct range with local measurement noise
        function p = p_hit(obj, z, z_deterministic)            
            if (z >= 0 && z <= obj.z_max && z_deterministic < obj.z_max)
                eta = normcdf(obj.z_max, z, obj.var_hit) - normcdf(0, z, obj.var_hit); % normalization factor (equation 6.6)
                p = normpdf(z, z_deterministic, obj.var_hit) / eta;
            else
                p = 0;
            end
        end
        
        % Unexpected objects
        function p = p_short(obj, z, z_deterministic)
            if (z >= 0 && z <= z_deterministic && z <= obj.z_max)
                eta = 1 - exp(-obj.lambda_short * z_deterministic); % normalization factor (equation 6.8)
                p = obj.lambda_short * exp(-obj.lambda_short * z) / eta;
            else
                p = 0;
            end
        end
        
        % Failures
        function p = p_max(obj, z, z_deterministic)
            if (z == obj.z_max)
                p = 1;
            else
                p = 0;
            end
        end       
        
        % Random measurements
        function p = p_rand(obj, z, z_deterministic)
            if (z >= 0 && z <= obj.z_max)
                p = 1 / obj.z_max;
            else
                p = 0;
            end
        end       
        
    end
end