classdef Dubins
    % See http://planning.cs.uiuc.edu/node821.html

    % Other repositories:
    % https://github.com/EwingKang/Dubins-Curve-For-MATLAB
    % https://github.com/AndrewWalker/Dubins-Curves#shkel01
    
    % OBS. Most of the code is flexible and supports different curvatures
    % although the original/traditional Dubins path consists of maximum
    % curvature and/or straight line segments - not variable curvature!

    % See also https://github.com/FelicienC/RRT-Dubins/blob/master/code/dubins.py
    
    properties (SetAccess = protected)
        q0 % initial configuration space (SE2 = [x,y,theta])
        q1 % terminal configuration space (SE2 = [x,y,theta])
        word % dubins transition word
        length % distance (total length) 
        q % intermediate poses
        params
        lengths
    end
    
    properties (Access = protected)                
    end
    
    methods (Static)
        function obj = from_word(q0, word, params)
            obj = Dubins();
            obj.q0 = q0;
            
            % Valid words:
            % LRL, RLR, LSL, LSR, RSL, RSR
            if (length(word) ~= 3)
                error('Incorrect word');
            end
            
            if (~contains('RL', word(1)))
                error('Incorrect word');
            end
            if (~contains('RLS', word(2)))
                error('Incorrect word');
            end
            if (~contains('RL', word(3)))
                error('Incorrect word');
            end            
            if (word(1) == word(2))
                error('Incorrect word');
            end
            if (word(2) == word(3))
                error('Incorrect word');
            end  
            
            obj.word = word;
            
            obj.params = params;
            obj = obj.compute_lengths();
        end 
    end
    
    methods                
        function obj = Dubins()            
        end       
        
        function draw(obj, varargin)
            if (~isempty(varargin))
                q1 = Dubins.draw_turn(obj.q0, obj.params(1), obj.params(2), varargin{1});
                if (obj.word(2) == 'S')
                    q2 = Dubins.draw_straight(q1, obj.params(3), varargin{1});
                else
                    q2 = Dubins.draw_turn(q1, obj.params(3), obj.params(4), varargin{1});
                end
                q3 = Dubins.draw_turn(q2, obj.params(5), obj.params(6), varargin{1});
            else
                q1 = Dubins.draw_turn(obj.q0, obj.params(1), obj.params(2));
                if (obj.word(2) == 'S')
                    q2 = Dubins.draw_straight(q1, obj.params(3));
                else
                    q2 = Dubins.draw_turn(q1, obj.params(3), obj.params(4));
                end
                q3 = Dubins.draw_turn(q2, obj.params(5), obj.params(6));
            end
        end          
        
        function obj = update_initial(obj, q0)
            obj.q0 = q0;
            obj.q1 = com
        end
        
        function q = get_pose(obj, s)
            if (any(s < 0) || any(s > obj.length))
                error('Incorrect query distance');
            end
            
            q = zeros(3, length(s));                        
            for (i = 1:length(s))
                if (s(i) < obj.lengths(1))
                    progress = s(i) / obj.lengths(1);
                    q(:,i) = Dubins.get_pose_turn(obj.q0, obj.params(1), obj.params(2) * progress);
                elseif (s(i) < (obj.lengths(1)+obj.lengths(2)))
                    progress = (s(i) - obj.lengths(1)) / obj.lengths(2);
                    if (obj.word(2) == 'S')
                        q(:,i) = Dubins.get_pose_straight(obj.q(:,2), obj.params(3) * progress);
                    else
                        q(:,i) = Dubins.get_pose_turn(obj.q(:,2), obj.params(3), obj.params(4) * progress);
                    end
                else
                    progress = (s(i) - obj.lengths(1) - obj.lengths(2)) / obj.lengths(3);
                    q(:,i) = Dubins.get_pose_turn(obj.q(:,3), obj.params(5), obj.params(6) * progress);
                end       
            end
        end 
    end
    
    methods (Access = private)           
        function obj = compute_lengths(obj)
            obj.lengths = zeros(3,1);
            
            obj.lengths(1) = Dubins.arc_length(obj.params(1), obj.params(2));
            
            if (obj.word(2) == 'S')
                obj.lengths(2) = obj.params(3);
            else
                obj.lengths(2) = Dubins.arc_length(obj.params(3), obj.params(4));
            end
            
            obj.lengths(3) = Dubins.arc_length(obj.params(5), obj.params(6));
            
            obj.length = sum(obj.lengths);            
            
            obj.q = zeros(3,4);
            obj.q(:,1) = obj.q0;
            obj.q(:,2) = Dubins.get_pose_turn(obj.q0, obj.params(1), obj.params(2));
            if (obj.word(2) == 'S')
                obj.q(:,3) = Dubins.get_pose_straight(obj.q(:,2), obj.params(3));
            else
                obj.q(:,3) = Dubins.get_pose_turn(obj.q(:,2), obj.params(3), obj.params(4));
            end
            obj.q(:,4) = Dubins.get_pose_turn(obj.q(:,3), obj.params(5), obj.params(6));            
            
            obj.q1 = obj.q(:,4);
        end
    end
    
    methods (Static)
        function obj = fit(q0, q1, varargin)
            if (~isempty(varargin))
                max_curvature = varargin{1};
            else
                max_curvature = 0.5;
            end
            
            [p1, length1] = Dubins.find_parameters_CSC(q0, q1, max_curvature);
            [p2, length2] = Dubins.find_parameters_CCC(q0, q1, max_curvature, false); % RLR
            [p3, length3] = Dubins.find_parameters_CCC(q0, q1, max_curvature, true);  % LRL

            if (length1 < length2 && length1 < length3)                
                params = [p1(1:3), 0, p1(4:5)];                              
            else        
                if (length2 < length3)
                    params = p2;
                else
                    params = p3;
                end                                
            end
            
            if (isempty(params))
                error('No solution was found!');
            end
            
            word = 'xxx';
            if (params(2) > 0)
                word(1) = 'L';
            elseif (params(2) < 0)
                word(1) = 'R';
            else
                word(1) = 'S';
            end
            if (params(4) > 0)
                word(2) = 'L';
            elseif (params(4) < 0)
                word(2) = 'R';
            else
                word(2) = 'S';
            end    
            if (params(6) > 0)
                word(3) = 'L';
            elseif (params(6) < 0)
                word(3) = 'R';
            else
                word(3) = 'S';
            end  
            
            obj = Dubins.from_word(q0, word, params);
            obj.q1 = q1;
        end  
        
        function obj = fit2(q0, q1, varargin)
            if (~isempty(varargin))
                max_curvature = varargin{1};
            else
                max_curvature = 0.5;
            end
            
            % When the distance between the points in greater than four 
            % times the minimum turn radius, the shortest path could be one
            % of the four combinations of CSC paths 
            r = 1/max_curvature;
            delta = q1(1:2) - q0(1:2);
            if (norm(delta) > 4*r)    
                % Find solution among CSC paths                
                lsr = Dubins.find_lsr(q0, q1, max_curvature);
                rsl = Dubins.find_rsl(q0, q1, max_curvature);
                lsl = Dubins.find_lsl(q0, q1, max_curvature);
                rsr = Dubins.find_rsr(q0, q1, max_curvature);

                paths = {lsr, rsl, lsl, rsr};
                lengths = [lsr.length, rsl.length, lsl.length, rsr.length];
            else
                % Find solution among CCC paths (RLR or LRL)
                rlr = Dubins.find_rlr(q0, q1, max_curvature);
                lrl = Dubins.find_lrl(q0, q1, max_curvature);
                
                paths = {rlr, lrl};
                lengths = [rlr.length, lrl.length];
            end
            
            [~,idx] = min(lengths);
            obj = paths{idx};
        end          
        
        function draw_csc(q0, c1,a1, d, c2,a2, varargin)
            if (~isempty(varargin))
                q1 = Dubins.draw_turn(q0, c1, a1, varargin{1});
                q2 = Dubins.draw_straight(q1, d, varargin{1});
                q3 = Dubins.draw_turn(q2, c2, a2, varargin{1});
            else
                q1 = Dubins.draw_turn(q0, c1, a1);
                q2 = Dubins.draw_straight(q1, d);
                q3 = Dubins.draw_turn(q2, c2, a2);
            end
        end
        
        function draw_ccc(q0, c1,a1, c2,a2, c3,a3, varargin)
            if (~isempty(varargin))
                q1 = Dubins.draw_turn(q0, c1, a1, varargin{1});                
                q2 = Dubins.draw_turn(q1, c2, a2, varargin{1});
                q3 = Dubins.draw_turn(q2, c3, a3, varargin{1});
            else
                q1 = Dubins.draw_turn(q0, c1, a1);                
                q2 = Dubins.draw_turn(q1, c2, a2);
                q3 = Dubins.draw_turn(q2, c3, a3);
            end
        end        
    end
    
    methods (Static)%, Access = private)
        function q = draw_straight(q, distance, varargin)
            if (~isempty(varargin))
                color = varargin{1};
            else
                color = 'b';
            end
            
            x0 = q(1);
            y0 = q(2);
            theta0 = q(3);
                        
            % Straight line
            x1 = x0 + distance * cos(theta0);
            y1 = y0 + distance * sin(theta0);
            line([x0;x1], [y0;y1], 'Color', color);
            
            q = [x1;y1;theta0];
        end
        
        function q = draw_turn(q, curvature, alpha, varargin)
            if (~isempty(varargin))
                color = varargin{1};
            else
                color = 'b';
            end
            
            x0 = q(1);
            y0 = q(2);
            theta0 = q(3);
            theta1 = theta0 + alpha; % ending angle
            
            % Arc radius
            R = 1/curvature;           
            
            % Depending on turning direction the arc center is on either
            % side of the trajectory
            if (alpha > 0)
                offset = -pi/2;                
            elseif (alpha < 0)
                offset = pi/2;
            else                
                return;
            end
            
            % Compute arc center
            xc = x0 - R*cos(theta0 + offset);
            yc = y0 - R*sin(theta0 + offset);
            
            arc_angles = 0:deg2rad(sign(alpha)*0.1):alpha; 

            arc_x = xc + R*cos(theta0 + arc_angles + offset);
            arc_y = yc + R*sin(theta0 + arc_angles + offset);

            line(arc_x, arc_y, 'Color', color);   
            
            q = [arc_x(end); arc_y(end); theta1];
        end  
        
        function q = get_pose_straight(q, d)
            x0 = q(1);
            y0 = q(2);
            theta0 = q(3);            

            x1 = x0 + d * cos(theta0);
            y1 = y0 + d * sin(theta0);            
            
            q = [x1;y1;theta0];
        end        
        
        function q = get_pose_turn(q, c, a)
            if (c < 0)
                error('Curvature can not be negative');
            end
            x0 = q(1);
            y0 = q(2);
            theta0 = q(3);            

            x1 = x0 + 1/c * (cos(theta0 + a - sign(a)*pi/2) - cos(theta0 - sign(a)*pi/2));
            y1 = y0 + 1/c * (sin(theta0 + a - sign(a)*pi/2) - sin(theta0 - sign(a)*pi/2));
            theta1 = theta0 + a; % ending angle
            
            q = [x1;y1;theta1];
        end
        
        function q = get_pose_turn2(q, c, a)
            if (c < 0)
                error('Curvature can not be negative');
            end
            x0 = q(1);
            y0 = q(2);
            theta0 = q(3);            

            x1 = x0 + 1/c * sign(a) * ( sin(theta0 + a) - sin(theta0));
            y1 = y0 + 1/c * sign(a) * (-cos(theta0 + a) + cos(theta0));
            theta1 = theta0 + a; % ending angle
            
            q = [x1;y1;theta1];
        end        
        
        function [p, dubins_length] = find_parameters_CSC(q0, q1, max_curvature)
            % Turn, straight, turn
            final_x = q1(1);
            final_y = q1(2);
            theta_delta = q1(3) - q0(3);
            pos_delta = norm(q1(1:2) - q0(1:2));
            
            % alpha1 + alpha2 = theta_delta
            % Optimization problem            
            % min arc_length(c1,a1) + d + arc_length(c2,a2)
            % find c1,a1, d, c2,a2
            % s.t. 
            %      dubins_x(q0, c1,a1, d, c2,a2) = final_x
            %      dubins_y(q0, c1,a1, d, c2,a2) = final_y
            %      a1+a2 = theta_delta
            %      -2*pi < a1 < 2*pi
            %      -2*pi < a2 < 2*pi
            %      c1 > 0
            %      c2 > 0
            %      d > 0
            
            import casadi.*;
            x0 = q0(1);
            y0 = q0(2);
            theta0 = q0(3);
            
            c1 = casadi.SX.sym('c1');
            a1 = casadi.SX.sym('a1');
            d = casadi.SX.sym('d');
            c2 = casadi.SX.sym('c2');
            a2 = casadi.SX.sym('a2');
            
            [x1,y1,theta1] = Dubins.dubins_csc_compact([x0;y0;theta0], c1,a1, d, c2,a2);
            
            nlp = struct('x', [c1,a1, d, c2,a2]', ...
                         'f', Dubins.arc_length(c1,a1) + d + Dubins.arc_length(c2,a2), ...
                         'g', [x1 - final_x;
                               y1 - final_y;
                               a1+a2 - theta_delta;
                               a1;
                               a2;
                               c1;
                               c2;
                               d]);
            S = casadi.nlpsol('S', 'ipopt', nlp);
            
            lb = [0;
                  0;
                  0;
                  -2*pi;
                  -2*pi;
                  0;  % min. curvature
                  0;  % min. curvature
                  0];     % min. distance
              
            ub = [0;
                  0;
                  0;
                  2*pi;
                  2*pi;
                  max_curvature;  % max. curvature
                  max_curvature;  % max. curvature 
                  inf]; % max. distance
              
            guess = [1,theta_delta/2, pos_delta, 1,theta_delta/2];
            
           
            r = S('x0', guess, 'lbg', lb, 'ubg', ub);
            x_opt = r.x;
            disp(x_opt)
            
            if (S.stats().success)
                p = zeros(1, length(x_opt));
                for (i = 1:length(x_opt))
                    p(i) = to_double(x_opt(i));
                end

                dubins_length = to_double(r.f);
            else
                p = [];
                dubins_length = inf;
            end
        end
        
        function [p, dubins_length] = find_parameters_CCC(q0, q1, max_curvature, LRL)
            % Turn, turn, turn
            final_x = q1(1);
            final_y = q1(2);
            theta_delta = q1(3) - q0(3);
            
            % alpha1 + alpha2 + alpha3 = theta_delta
            % Optimization problem            
            % min arc_length(c1,a1) + arc_length(c2,a2) + arc_length(c3,a3)
            % find c1,a1, c2,a2, c3,a3
            % s.t. 
            %      dubins_x(q0, c1,a1, c2,a2, c3,a3) = final_x
            %      dubins_y(q0, c1,a1, c2,a2, c3,a3) = final_y
            %      a1+a2+a3 = theta_delta
            %      -2*pi < a1 < 2*pi
            %      -2*pi < a2 < 2*pi      % this should maybe rather have been pi < a2 < 2*pi or -2*pi < a2 < -pi
            %      -2*pi < a3 < 2*pi
            %      c1 > 0
            %      c2 > 0
            %      c3 > 0            
            
            import casadi.*;
            x0 = q0(1);
            y0 = q0(2);
            theta0 = q0(3);
            
            c1 = casadi.SX.sym('c1');
            a1 = casadi.SX.sym('a1');            
            c2 = casadi.SX.sym('c2');
            a2 = casadi.SX.sym('a2');
            c3 = casadi.SX.sym('c3');
            a3 = casadi.SX.sym('a3');            
            
            [x1,y1,theta1] = Dubins.dubins_ccc_compact([x0;y0;theta0], c1,a1, c2,a2, c3,a3);
            
            nlp = struct('x', [c1,a1, c2,a2, c3,a3]', ...
                         'f', Dubins.arc_length(c1,a1) + Dubins.arc_length(c2,a2) + Dubins.arc_length(c3,a3), ...
                         'g', [x1 - final_x;
                               y1 - final_y;
                               a1+a2+a3 - theta_delta;
                               a1;
                               a2;
                               a3;
                               c1;
                               c2;
                               c3]);
            S = casadi.nlpsol('S', 'ipopt', nlp);
            
            lb = [0;
                  0;
                  0;
                  -2*pi;
                  -2*pi;
                  -2*pi;
                  0;  % min. curvature
                  0;  % min. curvature
                  0]; % min. curvature                  
              
            ub = [0;
                  0;
                  0;
                  2*pi;
                  2*pi;
                  2*pi;
                  max_curvature;  % max. curvature
                  max_curvature;  % max. curvature 
                  max_curvature];  % max. curvature                   
              
            if (LRL)
                lb(4) = 0;                
                ub(5) = -pi;                                              
                lb(6) = 0;
            else
                ub(4) = 0;               
                lb(5) = pi;              
                ub(6) = 0;
            end
              
            guess = [max_curvature,theta_delta/2, max_curvature,0, max_curvature,theta_delta/2];
           
            r = S('x0', guess, 'lbg', lb, 'ubg', ub);
            x_opt = r.x;
            disp(x_opt)
            
            if (S.stats().success)
                p = zeros(1, length(x_opt));
                for (i = 1:length(x_opt))
                    p(i) = to_double(x_opt(i));
                end
    
                dubins_length = to_double(r.f);
            else
                p = [];
                dubins_length = inf;
            end
        end        
        
        function l = arc_length(c,a)
            %R = 1/c;
            %l = 2*pi*R * abs(a)/(2*pi);
            l = abs(a) / c;
        end
        
        function [x,y,theta] = dubins_csc(q0, c1,a1, d, c2,a2)
            x = q0(1);
            y = q0(2);
            theta = q0(3);
            
            %% First section
            % Depending on turning direction the arc center is on either
            % side of the trajectory
            if (a1 > 0)
                offset1 = -pi/2;                
            else
                offset1 = pi/2;
            end
            x = x + 1/c1 * (cos(theta + a1 + offset1) - cos(theta + offset1));
            y = y + 1/c1 * (sin(theta + a1 + offset1) - sin(theta + offset1));                        
            theta = theta + a1;
            
            %% Second section
            % Straight line
            x = x + d * cos(theta);
            y = y + d * sin(theta);
            
            %% Third section
            % Depending on turning direction the arc center is on either
            % side of the trajectory
            if (a2 > 0)
                offset2 = -pi/2;                
            else
                offset2 = pi/2;
            end
            x = x + 1/c2 * (cos(theta + a2 + offset2) - cos(theta + offset2));
            y = y + 1/c2 * (sin(theta + a2 + offset2) - sin(theta + offset2));                                    
            theta = theta + a2;        
        end
        
        function [x,y,theta] = dubins_csc_compact(q0, c1,a1, d, c2,a2)
            x0 = q0(1);
            y0 = q0(2);
            theta0 = q0(3);

            x = x0 ...
              + 1/c1 * (cos(theta0 + a1 - sign(a1)*pi/2) - cos(theta0 - sign(a1)*pi/2)) ...
              + d * cos(theta0 + a1) ...
              + 1/c2 * (cos(theta0 + a1 + a2 - sign(a2)*pi/2) - cos(theta0 + a1 - sign(a2)*pi/2));
          
            y = y0 ...
              + 1/c1 * (sin(theta0 + a1 - sign(a1)*pi/2) - sin(theta0 - sign(a1)*pi/2)) ...
              + d * sin(theta0 + a1) ...
              + 1/c2 * (sin(theta0 + a1 + a2 - sign(a2)*pi/2) - sin(theta0 + a1 - sign(a2)*pi/2));
          
            theta = theta0 + a1 + a2;            
        end     
        
        function [x,y,theta] = dubins_ccc_compact(q0, c1,a1, c2,a2, c3,a3)
            x0 = q0(1);
            y0 = q0(2);
            theta0 = q0(3);

            x = x0 ...
              + 1/c1 * (cos(theta0 + a1 - sign(a1)*pi/2) - cos(theta0 - sign(a1)*pi/2)) ...              
              + 1/c2 * (cos(theta0 + a1 + a2 - sign(a2)*pi/2) - cos(theta0 + a1 - sign(a2)*pi/2)) ...
              + 1/c3 * (cos(theta0 + a1 + a2 + a3 - sign(a3)*pi/2) - cos(theta0 + a1 + a2 - sign(a3)*pi/2));
          
            y = y0 ...
              + 1/c1 * (sin(theta0 + a1 - sign(a1)*pi/2) - sin(theta0 - sign(a1)*pi/2)) ...              
              + 1/c2 * (sin(theta0 + a1 + a2 - sign(a2)*pi/2) - sin(theta0 + a1 - sign(a2)*pi/2)) ...
              + 1/c3 * (sin(theta0 + a1 + a2 + a3 - sign(a3)*pi/2) - sin(theta0 + a1 + a2 - sign(a3)*pi/2));
          
            theta = theta0 + a1 + a2 + a3;            
        end          
        
        %% Dubins closed form solutions
        % Based on "Dubins - On Curves of Minimal Length with a Constraint
        % on Average Curvature, and with Prescribed Initial and Terminal Positions and Tangents"
        % and https://arxiv.org/pdf/1804.07238.pdf
        %
        % OBS: When the distance between the start and goal points 
        % is greater than four times the minimum turn radius, the shortest 
        % path could be one of the four combinations of CSC paths 
        function path = find_lsl(q0, q1, max_curvature)
            r = 1/max_curvature;
            
            x0 = q0(1);
            y0 = q0(2);
            theta0 = q0(3);
            
            x1 = q1(1);
            y1 = q1(2);
            theta1 = q1(3);
            
            dx = x1 - x0;
            dy = y1 - y0;
            dtheta = theta1 - theta0;
            
            % Optimal parameters leading to shortest path
            Ls = sqrt( (dx - r*sin(dtheta))^2 + (dy + r*cos(dtheta) - r)^2 );            
            phi1 = mod(atan2( (dy + r*cos(dtheta) - r), (dx - r*sin(dtheta)) ), 2*pi);
            phi2 = mod(dtheta - phi1, 2*pi);
            
            % Construct Dubins object with these parameters
            path = Dubins.from_word(q0, 'LSL', [max_curvature, phi1, Ls, 0, max_curvature, phi2]);
        end        

        function path = find_rsr(q0, q1, max_curvature)
            r = 1/max_curvature;
            
            x0 = q0(1);
            y0 = q0(2);
            theta0 = q0(3);
            
            x1 = q1(1);
            y1 = q1(2);
            theta1 = q1(3);
            
            dx = x1 - x0;
            dy = -(y1 - y0);
            dtheta = -(theta1 - theta0);
            
            % Optimal parameters leading to shortest path
            Ls = sqrt( (dx - r*sin(dtheta))^2 + (dy + r*cos(dtheta) - r)^2 );            
            phi1 = mod(atan2( (dy + r*cos(dtheta) - r), (dx - r*sin(dtheta)) ), 2*pi);
            phi2 = mod(dtheta - phi1, 2*pi);
            
            % Construct Dubins object with these parameters
            path = Dubins.from_word(q0, 'RSR', [max_curvature, -phi1, Ls, 0, max_curvature, -phi2]);
        end
        
        function path = find_rsl(q0, q1, max_curvature)
            r = 1/max_curvature;
            
            x0 = q0(1);
            y0 = q0(2);
            theta0 = q0(3);
            
            x1 = q1(1);
            y1 = q1(2);
            theta1 = q1(3);
            
            dx = -(y1 - y0);
            dy = x1 - x0;
            dtheta = theta1 - theta0 + pi/2;
            
            % Optimal parameters leading to shortest path
            Lcc = sqrt( (dx - r*sin(dtheta) - r)^2 + (dy + r*cos(dtheta))^2 );            
            if (Lcc^2 < 4*r^2)
                warning('Can not find a valid RSL transition');
                path = Dubins();
                return;
            end            
            Ls = sqrt( Lcc^2 - 4*r^2 );
            
            psi1 = atan2( dy + r*cos(dtheta), dx - r*sin(dtheta) - r );
            psi2 = atan2( 2*r, Ls );
            
            phi1 = mod(-psi1 + psi2 + pi/2, 2*pi);
            phi2 = mod(dtheta + phi1 - pi/2, 2*pi);
            
            % Construct Dubins object with these parameters
            path = Dubins.from_word(q0, 'RSL', [max_curvature, -phi1, Ls, 0, max_curvature, phi2]);
        end        
        
        function path = find_lsr(q0, q1, max_curvature)
            r = 1/max_curvature;
            
            x0 = q0(1);
            y0 = q0(2);
            theta0 = q0(3);
            
            x1 = q1(1);
            y1 = q1(2);
            theta1 = q1(3);
            
            dx = (y1 - y0);
            dy = (x1 - x0);
            dtheta = -(theta1 - theta0) + pi/2;
            
            % Optimal parameters leading to shortest path
            Lcc = sqrt( (dx - r*sin(dtheta) - r)^2 + (dy + r*cos(dtheta))^2 );            
            if (Lcc^2 < 4*r^2)
                warning('Can not find a valid LSR transition');
                path = Dubins();
                return;
            end
            Ls = sqrt( Lcc^2 - 4*r^2 );                        
            
            psi1 = atan2( dy + r*cos(dtheta), dx - r*sin(dtheta) - r );
            psi2 = atan2( 2*r, Ls );
            
            phi1 = mod(-psi1 + psi2 + pi/2, 2*pi);
            phi2 = mod(dtheta + phi1 - pi/2, 2*pi);
            
            % Construct Dubins object with these parameters
            path = Dubins.from_word(q0, 'LSR', [max_curvature, phi1, Ls, 0, max_curvature, -phi2]);
        end      
        
        function path = find_lrl(q0, q1, max_curvature, varargin)
            if (~isempty(varargin))
                mid_circle_top = varargin{1}; % true or false
            else
                mid_circle_top = true;
            end
            r = 1/max_curvature;
            
            x0 = q0(1);
            y0 = q0(2);
            theta0 = q0(3);
            
            x1 = q1(1);
            y1 = q1(2);
            theta1 = q1(3);
            
%             dx = x1 - x0;
%             dy = y1 - y0;
%             dtheta = theta1 - theta0;
%             
%             Lcc = sqrt( dx^2 + dy^2 );    
%             mid = [dx;dy]/2;
%             
%             if (Lcc < 2*r)
%                 warning('Can not find a valid LRL transition');
%                 path = Dubins();
%                 return;
%             end
%            
%             psia = atan2(mid(2), mid(1));
%             
%             gamma = 2*asin(Lcc/(4*r));
%             phi1 = mod( (psia-theta0+pi/2+(pi-gamma)/2), (2*pi) );
%             phi3 = mod( (-psia+pi/2+theta1+(pi-gamma)/2), (2*pi) );
%             phi2 = 2*pi - gamma;
%             %total_len = (2*np.pi-gamma+abs(beta_0)+abs(beta_1))*r;
            
            c1 = q0(1:2)' + r*[-sin(q0(3)); cos(q0(3))];
            c3 = q1(1:2)' + r*[-sin(q1(3)); cos(q1(3))];

            v = c3 - c1;
            D = norm(v);
            if (D > 4*r)
                warning('Can not find a valid LRL transition');
                path = Dubins();
                return;
            end            
            % absolute angle to the x axis of the vector between first
            % and second turn centers
            if (mid_circle_top)
                gamma = atan2(v(2), v(1)) - acos(D / (4*r));
            else
                gamma = atan2(v(2), v(1)) + acos(D / (4*r));
            end

            c2 = c1 + 2*r*[cos(gamma); sin(gamma)];
            u = c2 - c3;
            psi = atan2(u(2), u(1));

            phi1 = gamma + pi/2 - theta0;
            phi2 = gamma - psi;
            phi3 = 3*pi/2 + theta1 - psi;
            
            phi1 = mod(phi1, 2*pi);
            phi2 = mod(phi2, 2*pi);
            phi3 = mod(phi3, 2*pi);
            
            path = Dubins.from_word(q0, 'LRL', [max_curvature, phi1, max_curvature, -phi2, max_curvature, phi3]);
        end
        
        function path = find_rlr(q0, q1, max_curvature, varargin)
            if (~isempty(varargin))
                mid_circle_top = varargin{1}; % true or false
            else
                mid_circle_top = true;
            end
            r = 1/max_curvature;
            
            x0 = q0(1);
            y0 = q0(2);
            theta0 = q0(3);
            
            x1 = q1(1);
            y1 = q1(2);
            theta1 = q1(3);
            
            dx = x1 - x0;
            dy = y1 - y0;
            dtheta = theta1 - theta0;
            
%             Lcc = sqrt( dx^2 + dy^2 );    
%             mid = [dx;dy]/2;
%             
%             if (Lcc < 2*r)
%                 warning('Can not find a valid RLR transition');
%                 path = Dubins();
%                 return;
%             end
%             
%             psia = atan2(mid(2), mid(1));
%             
%             gamma = 2*asin(Lcc/(4*r));
%             phi1 = -mod( (-psia+(theta0+pi/2)+(pi-gamma)/2), (2*pi) );
%             phi3 = -mod( (psia+pi/2-theta1+(pi-gamma)/2), (2*pi) );
%             phi2 = 2*pi - gamma;
%             %total_len = (2*np.pi-gamma+abs(beta_0)+abs(beta_1))*r;
%            

            c1 = q0(1:2)' - r*[-sin(q0(3)); cos(q0(3))];
            c3 = q1(1:2)' - r*[-sin(q1(3)); cos(q1(3))];

            v = c3 - c1;
            D = norm(v);
            if (D > 4*r)
                warning('Can not find a valid LRL transition');
                path = Dubins();
                return;
            end            
            % absolute angle to the x axis of the vector between first
            % and second turn centers            
            if (mid_circle_top)
                gamma = atan2(v(2), v(1)) - acos(D / (4*r));
            else
                gamma = atan2(v(2), v(1)) + acos(D / (4*r));
            end

            c2 = c1 + 2*r*[cos(gamma); sin(gamma)];
            u = c2 - c3;
            psi = atan2(u(2), u(1));

            phi1 = theta0 + pi/2 - gamma;
            phi2 = psi - gamma;
            phi3 = 3*pi/2 + psi - theta1;
            
            phi1 = mod(phi1, 2*pi);
            phi2 = mod(phi2, 2*pi);
            phi3 = mod(phi3, 2*pi);
            
            path = Dubins.from_word(q0, 'RLR', [max_curvature, -phi1, max_curvature, phi2, max_curvature, -phi3]);
        end        
        
    end
end