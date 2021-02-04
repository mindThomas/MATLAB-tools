classdef Robot
    properties (SetAccess = private)
        length
        width        
        pose % [x, y, heading]
    end
    
    methods
        function obj = Robot(length, width, init_pose)            
            obj.length = length;
            obj.width = width;
            obj.pose = init_pose;
        end
        
        function draw(obj)
            Robot.draw_(obj.pose(1), obj.pose(2), obj.pose(3), obj.length, obj.width);
        end
    end
               
    methods (Static)
        function c = get_corners(x, y, heading, length, width)
            car_corners = [
                length/2, -width/2; % front right
                length/2, width/2;  % front left            
                -length/2, width/2;  % rear left
                -length/2, -width/2; % rear right
            ];

            R = [cos(heading), -sin(heading);
                 sin(heading),  cos(heading)];
            
            oriented_corners = (R * car_corners')'; 
            
            c = [x, y] + oriented_corners;   
        end
        
        function draw_(x, y, heading, length, width, varargin)             
            if (~isempty(varargin))
                color = varargin{1};
            else
                color = 'b';
            end
            
            c = Robot.get_corners(x, y, heading, length, width);
            plot([c(1,1),c(2,1),c(3,1),c(4,1),c(1,1)], ...
                            [c(1,2),c(2,2),c(3,2),c(4,2),c(1,2)], color);
        end          
        
        function pose_next = step(pose, c, d)
            % Drive in a constant curvature, c, with arc curve length, d   
            x0 = pose(1);
            y0 = pose(2);
            theta0 = pose(3);
            %if (abs(c) > 0.001)
                % convert arc curve length into angle
                a = d .* c;
                %x1 = x0 + 1./abs(c) .* (cos(theta0 + a - sign(a)*pi/2) - cos(theta0 - sign(a)*pi/2));
                %y1 = y0 + 1./abs(c)  .* (sin(theta0 + a - sign(a)*pi/2) - sin(theta0 - sign(a)*pi/2));
                x1 = x0 + 1./abs(c) .* sign(a) .* (sin(theta0 + a) - sin(theta0));
                y1 = y0 + 1./abs(c)  .* sign(a) .* (-cos(theta0 + a) + cos(theta0));
                theta1 = theta0 + a; % ending angle
            %else % drive straight
                x1(abs(c) <= 0.001) = x0 + d(abs(c) <= 0.001) * cos(theta0);
                y1(abs(c) <= 0.001) = y0 + d(abs(c) <= 0.001) * sin(theta0);                    
                theta1(abs(c) <= 0.001) = theta0;
            %end
            
            pose_next = zeros(length(c), 3);
            pose_next(:,1) = x1;
            pose_next(:,2) = y1;
            pose_next(:,3) = theta1;
        end
        
        function collision = in_collision(x, y, heading, length, width, grid, ds)
            c = Robot.get_corners(x, y, heading, length, width);
                        
            for (i = 1:size(c, 1))
                % Assuming bottom left corner of grid is (0,0)
                ix = floor(c(i,1) / ds) + 1;
                iy = size(grid, 1) - floor(c(i,2) / ds);
                
                for (dx = 0:1)
                    for (dy = 0:1)
                        iix = ix + dx;
                        iiy = iy + dy;
                        if (iix > 0 && iix <= size(grid, 2) && ...
                            iiy > 0 && iiy <= size(grid, 1))
                            if (~grid(iiy, iix))
                                collision = true;
                                return;
                            end
                        end
                    end
                end
            end
            
            collision = false;
        end
    end
end