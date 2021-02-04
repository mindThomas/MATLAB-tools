classdef ReedsShepp
    % See http://planning.cs.uiuc.edu/node822.html
    %
    % The Reeds-Shepp curves are extensions of Dubins curves.
    % Dubins curves are applied when finding the shortest trajectory 
    % between two points in R^2. If there is no limitations on the 
    % steering angle, then the shortest path between start and goal is 
    % a straight line. In Dubins curve the speed is constant and the 
    % vehicle only has a forward direction.
    %
    % Based on: https://github.com/ghliu/pyReedsShepp/blob/master/reeds_shepp/src/reeds_shepp.cpp
    
    properties (SetAccess = protected)
        q0 % initial configuration space (SE2 = [x,y,theta])
        q1 % terminal configuration space (SE2 = [x,y,theta])
        type % Reeds Shepp transition word        
        length % distance (in Reeds Shepp unit frame)
        q % intermediate poses        
        lengths
        params
        
        rho % turning radius
    end
    
    properties (Constant)
        RS_NOP=' ';
        RS_LEFT='L';
        RS_STRAIGHT='S';
        RS_RIGHT='R';
     
        PATH_TYPES = [...
            [ReedsShepp.RS_LEFT, ReedsShepp.RS_RIGHT, ReedsShepp.RS_LEFT, ReedsShepp.RS_NOP, ReedsShepp.RS_NOP]             % 0
            [ReedsShepp.RS_RIGHT, ReedsShepp.RS_LEFT, ReedsShepp.RS_RIGHT, ReedsShepp.RS_NOP, ReedsShepp.RS_NOP]            % 1
            [ReedsShepp.RS_LEFT, ReedsShepp.RS_RIGHT, ReedsShepp.RS_LEFT, ReedsShepp.RS_RIGHT, ReedsShepp.RS_NOP]           % 2
            [ReedsShepp.RS_RIGHT, ReedsShepp.RS_LEFT, ReedsShepp.RS_RIGHT, ReedsShepp.RS_LEFT, ReedsShepp.RS_NOP]           % 3
            [ReedsShepp.RS_LEFT, ReedsShepp.RS_RIGHT, ReedsShepp.RS_STRAIGHT, ReedsShepp.RS_LEFT, ReedsShepp.RS_NOP]        % 4
            [ReedsShepp.RS_RIGHT, ReedsShepp.RS_LEFT, ReedsShepp.RS_STRAIGHT, ReedsShepp.RS_RIGHT, ReedsShepp.RS_NOP]       % 5
            [ReedsShepp.RS_LEFT, ReedsShepp.RS_STRAIGHT, ReedsShepp.RS_RIGHT, ReedsShepp.RS_LEFT, ReedsShepp.RS_NOP]        % 6
            [ReedsShepp.RS_RIGHT, ReedsShepp.RS_STRAIGHT, ReedsShepp.RS_LEFT, ReedsShepp.RS_RIGHT, ReedsShepp.RS_NOP]       % 7
            [ReedsShepp.RS_LEFT, ReedsShepp.RS_RIGHT, ReedsShepp.RS_STRAIGHT, ReedsShepp.RS_RIGHT, ReedsShepp.RS_NOP]       % 8
            [ReedsShepp.RS_RIGHT, ReedsShepp.RS_LEFT, ReedsShepp.RS_STRAIGHT, ReedsShepp.RS_LEFT, ReedsShepp.RS_NOP]        % 9
            [ReedsShepp.RS_RIGHT, ReedsShepp.RS_STRAIGHT, ReedsShepp.RS_RIGHT, ReedsShepp.RS_LEFT, ReedsShepp.RS_NOP]       % 10
            [ReedsShepp.RS_LEFT, ReedsShepp.RS_STRAIGHT, ReedsShepp.RS_LEFT, ReedsShepp.RS_RIGHT, ReedsShepp.RS_NOP]        % 11
            [ReedsShepp.RS_LEFT, ReedsShepp.RS_STRAIGHT, ReedsShepp.RS_RIGHT, ReedsShepp.RS_NOP, ReedsShepp.RS_NOP]         % 12
            [ReedsShepp.RS_RIGHT, ReedsShepp.RS_STRAIGHT, ReedsShepp.RS_LEFT, ReedsShepp.RS_NOP, ReedsShepp.RS_NOP]         % 13
            [ReedsShepp.RS_LEFT, ReedsShepp.RS_STRAIGHT, ReedsShepp.RS_LEFT, ReedsShepp.RS_NOP, ReedsShepp.RS_NOP]          % 14
            [ReedsShepp.RS_RIGHT, ReedsShepp.RS_STRAIGHT, ReedsShepp.RS_RIGHT, ReedsShepp.RS_NOP, ReedsShepp.RS_NOP]        % 15
            [ReedsShepp.RS_LEFT, ReedsShepp.RS_RIGHT, ReedsShepp.RS_STRAIGHT, ReedsShepp.RS_LEFT, ReedsShepp.RS_RIGHT]      % 16
            [ReedsShepp.RS_RIGHT, ReedsShepp.RS_LEFT, ReedsShepp.RS_STRAIGHT, ReedsShepp.RS_RIGHT, ReedsShepp.RS_LEFT]...       % 17
        ];
    
        RS_EPS = 1e-6;
        ZERO = 10*eps;
    end
    
    properties (Access = protected)                
    end
    
    methods (Static)
        function path = reedsShepp(x, y, phi, backwards_enabled)
            path = ReedsShepp();
            path = ReedsShepp.CSC(x, y, phi, path, backwards_enabled);
            path = ReedsShepp.CCC(x, y, phi, path, backwards_enabled);
            path = ReedsShepp.CCCC(x, y, phi, path, backwards_enabled);
            path = ReedsShepp.CCSC(x, y, phi, path, backwards_enabled);
            %path = ReedsShepp.CCSCC(x, y, phi, path, backwards_enabled); % This is not working correctly
            
            if (path.length == inf)
                error('No path was found');
            end
        end
        
        function path = fit(q0, q1, max_curvature, varargin)
            if (~isempty(varargin))
                backwards_enabled = varargin{1};
            else
                backwards_enabled = true;
            end
            
            rho = 1/max_curvature;
            dx = q1(1) - q0(1);
            dy = q1(2) - q0(2);
            dth = q1(3) - q0(3);
            c = cos(q0(3));
            s = sin(q0(3));
            x = c*dx + s*dy;
            y = -s*dx + c*dy;
            path = ReedsShepp.reedsShepp(x/rho, y/rho, dth, backwards_enabled );
            path.rho = rho;
            path.q0 = q0;
            path.q1 = q1;
            path.length = rho * path.length;
        end

        
        function obj = from_word2(q0, word, directions, params)
            obj = ReedsShepp();
            obj.q0 = q0;
            
            if (length(word) ~= length(directions))
                error('Invalid directions');
            end
            
            % Valid words:
            % See Figure 15.10 from http://planning.cs.uiuc.edu/node822.html
            if (length(word) == 3)
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
                if (word(2) == 'S')
                    if (max(directions) - min(directions))
                        % Directions should all be the same
                        error('Invalid directions');
                    end
                end                
            elseif (length(word) == 4)
                if (word(1) == 'S')
                    error('Incorrect word');
                end
                if (word(4) == 'S')
                    error('Incorrect word');
                end
                if (sum(abs(diff(directions))) > 2)                
                    error('Invalid directions');
                end                
            elseif (length(word) == 5)
                if (word(3) ~= 'S')
                    error('Incorrect word');
                end
                
                d = abs(diff(directions));
                if (sum(d) > 2)                
                    error('Invalid directions');
                end                                  
                if (d(1) ~= 1 || d(2) ~= 0 || d(3) ~= 0 || d(4) ~= 1)                
                    error('Invalid directions');
                end                  
            else
                error('Incorrect word');
            end
            
            obj.word = word;
            obj.directions = directions;
            
            obj.params = params;
            obj = obj.compute_lengths();
        end 
    end
    
    methods                
        function obj = ReedsShepp(varargin)                         
            type = repmat(ReedsShepp.RS_NOP, 1, 5);  
            t = inf;
            u = 0;
            v = 0;
            w = 0;
            x = 0;
            
            if (length(varargin) >= 2)
                type = varargin{1};
            end
            if (length(varargin) >= 2)
                t = varargin{2};
            end
            if (length(varargin) >= 3)
                u = varargin{3};
            end
            if (length(varargin) >= 4)
                v = varargin{4};
            end
            if (length(varargin) >= 5)
                w = varargin{5};
            end  
            if (length(varargin) >= 6)
                w = varargin{6};
            end              
                      
            obj.type = type;
            obj.params = zeros(5,1);
            obj.params(1) = t;
            obj.params(2) = u;
            obj.params(3) = v;
            obj.params(4) = w;
            obj.params(5) = x;
            obj.length = abs(t) + abs(u) + abs(v) + abs(w) + abs(x); % total length            
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
            if (s < 0.0)
                s = 0.0;
            end
            if (s > obj.length)
                s = obj.length;
            end
            s = s / obj.rho;

            q = [0, 0, obj.q0(3)];

            for (i = 1:5)
                if (s <= 0)
                    break;
                end

                if (obj.params(i) < 0)        
                    v = max(-s, obj.params(i));
                    s = s + v;        
                else        
                    v = min(s, obj.params(i));
                    s = s - v;
                end

                phi = q(3);

                if (obj.type(i) == ReedsShepp.RS_LEFT)
                    q(1) = q(1) + ( sin(phi+v) - sin(phi));
                    q(2) = q(2) + (-cos(phi+v) + cos(phi));
                    q(3) = phi + v;
                elseif (obj.type(i) == ReedsShepp.RS_RIGHT)
                    q(1) = q(1) + (-sin(phi-v) + sin(phi));
                    q(2) = q(2) + ( cos(phi-v) - cos(phi));
                    q(3) = phi - v;
                elseif (obj.type(i) == ReedsShepp.RS_STRAIGHT)
                    q(1) = q(1) + (v * cos(phi));
                    q(2) = q(2) + (v * sin(phi));
                end
            end

            q(1) = q(1) * obj.rho + obj.q0(1);
            q(2) = q(2) * obj.rho + obj.q0(2);
        end
    end   
    
    methods (Static)        
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
        
        %% Fitting
        function v = mod2pi(x)    
            v = mod(x, 2*pi);
            if (v < -pi)
                v = v + 2*pi;
            else
                if (v > pi)
                    v = v - 2*pi;
                end
            end        
        end

        function [r, theta] = polar(x, y)    
            r = sqrt(x*x + y*y);
            theta = atan2(y, x);
        end
        function [tau, omega] = tauOmega(u, v, xi, eta, phi)    
            delta = ReedsShepp.mod2pi(u-v);
            A = sin(u) - sin(delta);
            B = cos(u) - cos(delta) - 1.;
            t1 = atan2(eta*A - xi*B, xi*A + eta*B);
            t2 = 2. * (cos(delta) - cos(v) - cos(u)) + 3;
            if (t2<0)
                tau = ReedsShepp.mod2pi(t1+pi);
            else
                tau = ReedsShepp.mod2pi(t1);
            end
            omega = ReedsShepp.mod2pi(tau - u + v - phi);
        end

        % formula 8.1 in Reeds-Shepp paper
        function [valid, t,u,v] = LpSpLp(x, y, phi)    
            u=0; t=0; v=0;
            [u,t] = ReedsShepp.polar(x - sin(phi), y - 1. + cos(phi));
            if (t >= -ReedsShepp.ZERO)        
                v = ReedsShepp.mod2pi(phi - t);
                if (v >= -ReedsShepp.ZERO)        
                    assert(abs(u*cos(t) + sin(phi) - x) < ReedsShepp.RS_EPS);
                    assert(abs(u*sin(t) - cos(phi) + 1 - y) < ReedsShepp.RS_EPS);
                    assert(abs(ReedsShepp.mod2pi(t+v - phi)) < ReedsShepp.RS_EPS);
                    valid = true;
                    return;
                end
            end
            valid = false;
            return;
        end

        %formula 8.2    
        function [valid, t,u,v] = LpSpRp(x, y, phi)    
            u=0; t=0; v=0;
            [u1, t1] = ReedsShepp.polar(x + sin(phi), y - 1. - cos(phi));
            u1 = u1*u1;
            if (u1 >= 4.)            
                u = sqrt(u1 - 4.);
                theta = atan2(2., u);
                t = ReedsShepp.mod2pi(t1 + theta);
                v = ReedsShepp.mod2pi(t - phi);
                assert(abs(2*sin(t) + u*cos(t) - sin(phi) - x) < ReedsShepp.RS_EPS);
                assert(abs(-2*cos(t) + u*sin(t) + cos(phi) + 1 - y) < ReedsShepp.RS_EPS);
                assert(abs(ReedsShepp.mod2pi(t-v - phi)) < ReedsShepp.RS_EPS);
                valid = (t>=-ReedsShepp.ZERO && v>=-ReedsShepp.ZERO);
                return
            end
            valid = false;
            return;
        end

        function path = CSC(x, y, phi, path, backwards_enabled)    
            Lmin = path.length();

            [valid, t, u, v] = ReedsShepp.LpSpLp(x, y, phi);
            L = abs(t) + abs(u) + abs(v);
            positive = ((t>=0) && (u>=0) && (v>=0)) || backwards_enabled;
            if (valid && Lmin > L && positive)        
                path = ReedsShepp(ReedsShepp.PATH_TYPES(15,:), t, u, v);
                Lmin = L;
            end
            
            [valid, t, u, v] = ReedsShepp.LpSpLp(-x, y, -phi); % timeflip
            L = abs(t) + abs(u) + abs(v);
            positive = ((-t>=0) && (-u>=0) && (-v>=0)) || backwards_enabled;
            if (valid && Lmin > L && positive)      
                path = ReedsShepp(ReedsShepp.PATH_TYPES(15,:), -t, -u, -v);
                Lmin = L;
            end            

            [valid, t, u, v] = ReedsShepp.LpSpLp(x, -y, -phi); % reflect
            L = abs(t) + abs(u) + abs(v);
            positive = ((t>=0) && (u>=0) && (v>=0)) || backwards_enabled;
            if (valid && Lmin > L && positive)  
                path = ReedsShepp(ReedsShepp.PATH_TYPES(16,:), t, u, v);
                Lmin = L;
            end
            
            [valid, t, u, v] = ReedsShepp.LpSpLp(-x, -y, phi); % timeflip + reflect
            L = abs(t) + abs(u) + abs(v);
            positive = ((-t>=0) && (-u>=0) && (-v>=0)) || backwards_enabled;
            if (valid && Lmin > L && positive)   
                path = ReedsShepp(ReedsShepp.PATH_TYPES(16,:), -t, -u, -v);
                Lmin = L;
            end            

            [valid, t, u, v] = ReedsShepp.LpSpRp(x, y, phi);
            L = abs(t) + abs(u) + abs(v);
            positive = ((t>=0) && (u>=0) && (v>=0)) || backwards_enabled;
            if (valid && Lmin > L && positive)  
                path = ReedsShepp(ReedsShepp.PATH_TYPES(13,:), t, u, v);
                Lmin = L;
            end        
            
            [valid, t, u, v] = ReedsShepp.LpSpRp(-x, y, -phi); % timeflip
            L = abs(t) + abs(u) + abs(v);
            positive = ((-t>=0) && (-u>=0) && (-v>=0)) || backwards_enabled;
            if (valid && Lmin > L && positive)  
                path = ReedsShepp(ReedsShepp.PATH_TYPES(13,:), -t, -u, -v);
                Lmin = L;
            end             
            
            [valid, t, u, v] = ReedsShepp.LpSpRp(x, -y, -phi); % reflect
            L = abs(t) + abs(u) + abs(v);
            positive = ((t>=0) && (u>=0) && (v>=0)) || backwards_enabled;
            if (valid && Lmin > L && positive)   
                path = ReedsShepp(ReedsShepp.PATH_TYPES(14,:), t, u, v);
                Lmin = L;
            end 
            
            [valid, t, u, v] = ReedsShepp.LpSpRp(-x, -y, phi); % timeflip + reflect
            L = abs(t) + abs(u) + abs(v);
            positive = ((-t>=0) && (-u>=0) && (-v>=0)) || backwards_enabled;
            if (valid && Lmin > L && positive)  
                path = ReedsShepp(ReedsShepp.PATH_TYPES(14,:), -t, -u, -v);
                Lmin = L;            
            end
        end


        % formula 8.3 / 8.4  *** TYPO IN PAPER ***
        function [valid, t,u,v] = LpRmL(x, y, phi) 
            u=0; t=0; v=0;
            xi = x - sin(phi);
            eta = y - 1. + cos(phi);        
            [u1, theta] = ReedsShepp.polar(xi, eta);
            if (u1 <= 4.)        
                u = -2.*asin(.25 * u1);
                t = ReedsShepp.mod2pi(theta + .5 * u + pi);
                v = ReedsShepp.mod2pi(phi - t + u);
                assert(abs(2*(sin(t) - sin(t-u)) + sin(phi) - x) < ReedsShepp.RS_EPS);
                assert(abs(2*(-cos(t) + cos(t-u)) - cos(phi) + 1 - y) < ReedsShepp.RS_EPS);
                assert(abs(ReedsShepp.mod2pi(t-u+v - phi)) < ReedsShepp.RS_EPS);            
                valid = (t>=-ReedsShepp.ZERO && u<=ReedsShepp.ZERO);
                return            
            end
            valid = false;
            return;
        end

        function path = CCC(x, y, phi, path, backwards_enabled)    
            Lmin = path.length();

            [valid, t, u, v] = ReedsShepp.LpRmL(x, y, phi);
            L = abs(t) + abs(u) + abs(v);
            positive = ((t>=0) && (u>=0) && (v>=0)) || backwards_enabled;
            if (valid && Lmin > L && positive)          
                path = ReedsShepp(ReedsShepp.PATH_TYPES(1,:), t, u, v);
                Lmin = L;
            end
            
            [valid, t, u, v] = ReedsShepp.LpRmL(-x, y, -phi); % timeflip
            L = abs(t) + abs(u) + abs(v);
            positive = ((-t>=0) && (-u>=0) && (-v>=0)) || backwards_enabled;
            if (valid && Lmin > L && positive)   
                path = ReedsShepp(ReedsShepp.PATH_TYPES(1,:), -t, -u, -v);
                Lmin = L;            
            end

            [valid, t, u, v] = ReedsShepp.LpRmL(x, -y, -phi); % reflect
            L = abs(t) + abs(u) + abs(v);
            positive = ((t>=0) && (u>=0) && (v>=0)) || backwards_enabled;
            if (valid && Lmin > L && positive)       
                path = ReedsShepp(ReedsShepp.PATH_TYPES(1,:), t, u, v);
                Lmin = L;
            end     
            
            [valid, t, u, v] = ReedsShepp.LpRmL(-x, -y, phi); % timeflip + reflect
            L = abs(t) + abs(u) + abs(v);
            positive = ((-t>=0) && (-u>=0) && (-v>=0)) || backwards_enabled;
            if (valid && Lmin > L && positive)  
                path = ReedsShepp(ReedsShepp.PATH_TYPES(1,:), -t, -u, -v);
                Lmin = L;
            end                     

            % backwards
            xb = x*cos(phi) + y*sin(phi);
            yb = x*sin(phi) - y*cos(phi);

            [valid, t, u, v] = ReedsShepp.LpRmL(xb, yb, phi);
            L = abs(t) + abs(u) + abs(v);
            positive = ((t>=0) && (u>=0) && (v>=0)) || backwards_enabled;
            if (valid && Lmin > L && positive)  
                path = ReedsShepp(ReedsShepp.PATH_TYPES(1,:), v, u, t);
                Lmin = L;
            end  
            
            [valid, t, u, v] = ReedsShepp.LpRmL(-xb, yb, -phi); % timeflip
            L = abs(t) + abs(u) + abs(v);
            positive = ((-t>=0) && (-u>=0) && (-v>=0)) || backwards_enabled;
            if (valid && Lmin > L && positive)      
                path = ReedsShepp(ReedsShepp.PATH_TYPES(1,:), -v, -u, -t);
                Lmin = L;            
            end

            [valid, t, u, v] = ReedsShepp.LpRmL(xb, -yb, -phi); % reflect
            L = abs(t) + abs(u) + abs(v);
            positive = ((t>=0) && (u>=0) && (v>=0)) || backwards_enabled;
            if (valid && Lmin > L && positive)         
                path = ReedsShepp(ReedsShepp.PATH_TYPES(1,:), v, u, t);
                Lmin = L;
            end           
            
            [valid, t, u, v] = ReedsShepp.LpRmL(-xb, -yb, phi); % timeflip + reflect
            L = abs(t) + abs(u) + abs(v);
            positive = ((-t>=0) && (-u>=0) && (-v>=0)) || backwards_enabled;
            if (valid && Lmin > L && positive)  
                path = ReedsShepp(ReedsShepp.PATH_TYPES(1,:), -v, -u, -t);
                Lmin = L;                
            end
        end


        % formula 8.7
        function [valid, t,u,v] = LpRupLumRm(x, y, phi) 
            u=0; t=0; v=0;
            xi = x + sin(phi);
            eta = y - 1. - cos(phi);
            rho = .25 * (2. + sqrt(xi*xi + eta*eta));
            if (rho <= 1.)        
                u = acos(rho);
                [t, v] = ReedsShepp.tauOmega(u, -u, xi, eta, phi);
                assert(abs(2*(sin(t)-sin(t-u)+sin(t-2*u))-sin(phi) - x) < ReedsShepp.RS_EPS);
                assert(abs(2*(-cos(t)+cos(t-u)-cos(t-2*u))+cos(phi)+1 - y) < ReedsShepp.RS_EPS);
                assert(abs(ReedsShepp.mod2pi(t-2*u-v - phi)) < ReedsShepp.RS_EPS);            
                valid = (t>=-ReedsShepp.ZERO && v<=ReedsShepp.ZERO);
                return            
            end
            valid = false;
            return;
        end
        % formula 8.8
        function [valid, t,u,v] = LpRumLumRp(x, y, phi) 
            u=0; t=0; v=0;
            xi = x + sin(phi);
            eta = y - 1. - cos(phi);
            rho = (20. - xi*xi - eta*eta) / 16.;
            if (rho>=0 && rho<=1)        
                u = -acos(rho);
                if (u >= -.5 * pi)            
                    [t, v] = ReedsShepp.tauOmega(u, u, xi, eta, phi);
                    assert(abs(4*sin(t)-2*sin(t-u)-sin(phi) - x) < ReedsShepp.RS_EPS);
                    assert(abs(-4*cos(t)+2*cos(t-u)+cos(phi)+1 - y) < ReedsShepp.RS_EPS);
                    assert(abs(ReedsShepp.mod2pi(t-v - phi)) < ReedsShepp.RS_EPS);                
                    valid = (t>=-ReedsShepp.ZERO && v>=-ReedsShepp.ZERO);
                    return            
                end
            end
            valid = false;
            return;
        end

        function path = CCCC(x, y, phi, path, backwards_enabled)    
            Lmin = path.length();

            if (backwards_enabled)
                [valid, t, u, v] = ReedsShepp.LpRupLumRm(x, y, phi);
                L = abs(t) + 2.*abs(u) + abs(v);
                if (valid && Lmin > L)        
                    path = ReedsShepp(ReedsShepp.PATH_TYPES(3,:), t, u, -u, v);
                    Lmin = L;
                end

                [valid, t, u, v] = ReedsShepp.LpRupLumRm(-x, y, -phi); % timeflip
                L = abs(t) + 2.*abs(u) + abs(v);
                if (valid && Lmin > L)        
                    path = ReedsShepp(ReedsShepp.PATH_TYPES(3,:), -t, -u, u, -v);
                    Lmin = L;
                end 

                [valid, t, u, v] = ReedsShepp.LpRupLumRm(x, -y, -phi); % reflect
                L = abs(t) + 2.*abs(u) + abs(v);
                if (valid && Lmin > L)        
                    path = ReedsShepp(ReedsShepp.PATH_TYPES(4,:), t, u, -u, v);
                    Lmin = L;
                end         

                [valid, t, u, v] = ReedsShepp.LpRupLumRm(-x, -y, phi); % timeflip + reflect
                L = abs(t) + 2.*abs(u) + abs(v);
                if (valid && Lmin > L)        
                    path = ReedsShepp(ReedsShepp.PATH_TYPES(4,:), -t, -u, u, -v);
                    Lmin = L;
                end     
            end

            [valid, t, u, v] = ReedsShepp.LpRumLumRp(x, y, phi);
            L = abs(t) + 2.*abs(u) + abs(v);
            positive = ((t>=0) && (u>=0) && (v>=0)) || backwards_enabled;
            if (valid && Lmin > L && positive)        
                path = ReedsShepp(ReedsShepp.PATH_TYPES(3,:), t, u, u, v);
                Lmin = L;
            end     
            
            [valid, t, u, v] = ReedsShepp.LpRumLumRp(-x, y, -phi); % timeflip
            L = abs(t) + 2.*abs(u) + abs(v);
            positive = ((-t>=0) && (-u>=0) && (-v>=0)) || backwards_enabled;
            if (valid && Lmin > L && positive)        
                path = ReedsShepp(ReedsShepp.PATH_TYPES(3,:), -t, -u, -u, -v);
                Lmin = L;
            end                 

            [valid, t, u, v] = ReedsShepp.LpRumLumRp(x, -y, -phi); % reflect
            L = abs(t) + 2.*abs(u) + abs(v);
            positive = ((t>=0) && (u>=0) && (v>=0)) || backwards_enabled;
            if (valid && Lmin > L && positive)        
                path = ReedsShepp(ReedsShepp.PATH_TYPES(4,:), t, u, u, v);
                Lmin = L;
            end
            
            [valid, t, u, v] = ReedsShepp.LpRumLumRp(-x, -y, phi); % timeflip + reflect
            L = abs(t) + 2.*abs(u) + abs(v);
            positive = ((-t>=0) && (-u>=0) && (-v>=0)) || backwards_enabled;
            if (valid && Lmin > L && positive)             
                path = ReedsShepp(ReedsShepp.PATH_TYPES(4,:), -t, -u, -u, -v);
                Lmin = L;
            end
        end

        % formula 8.9
        function [valid, t,u,v] = LpRmSmLm(x, y, phi) 
            u=0; t=0; v=0;
            xi = x - sin(phi);
            eta = y - 1. + cos(phi);        
            [rho, theta] = ReedsShepp.polar(xi, eta);
            if (rho >= 2.)        
                r = sqrt(rho*rho - 4.);
                u = 2. - r;
                t = ReedsShepp.mod2pi(theta + atan2(r, -2.));
                v = ReedsShepp.mod2pi(phi - .5*pi - t);
                assert(abs(2*(sin(t)-cos(t))-u*sin(t)+sin(phi) - x) < ReedsShepp.RS_EPS);
                assert(abs(-2*(sin(t)+cos(t))+u*cos(t)-cos(phi)+1 - y) < ReedsShepp.RS_EPS);
                assert(abs(ReedsShepp.mod2pi(t+pi/2+v-phi)) < ReedsShepp.RS_EPS);            
                valid = (t>=-ReedsShepp.ZERO && u<=ReedsShepp.ZERO && v<=ReedsShepp.ZERO);
                return            
            end
            valid = false;
            return;
        end

        %formula 8.10
        function [valid, t,u,v] = LpRmSmRm(x, y, phi) 
            u=0; t=0; v=0;
            xi = x + sin(phi);
            eta = y - 1. - cos(phi);        
            [rho, theta] = ReedsShepp.polar(-eta, xi);
            if (rho >= 2.)        
                t = theta;
                u = 2. - rho;
                v = ReedsShepp.mod2pi(t + .5*pi - phi);
                assert(abs(2*sin(t)-cos(t-v)-u*sin(t) - x) < ReedsShepp.RS_EPS);
                assert(abs(-2*cos(t)-sin(t-v)+u*cos(t)+1 - y) < ReedsShepp.RS_EPS);
                assert(abs(ReedsShepp.mod2pi(t+pi/2-v-phi)) < ReedsShepp.RS_EPS);            
                valid = (t>=-ReedsShepp.ZERO && u<=ReedsShepp.ZERO && v<=ReedsShepp.ZERO);
                return            
            end
            valid = false;
            return;
        end

        function path = CCSC(x, y, phi, path, backwards_enabled)    
            Lmin = path.length() - .5*pi;

            if (~backwards_enabled)
                return; % all CCSC paths involves backwards driving
            end
            
            [valid, t, u, v] = ReedsShepp.LpRmSmLm(x, y, phi);
            L = abs(t) + abs(u) + abs(v);
            if (valid && Lmin > L)        
                path = ReedsShepp(ReedsShepp.PATH_TYPES(5,:), t, -.5*pi, u, v);
                Lmin = L;
            end

            [valid, t, u, v] = ReedsShepp.LpRmSmLm(-x, y, -phi); % timeflip
            L = abs(t) + abs(u) + abs(v);
            if (valid && Lmin > L)        
                path = ReedsShepp(ReedsShepp.PATH_TYPES(5,:), -t, .5*pi, -u, -v);
                Lmin = L;
            end

            [valid, t, u, v] = ReedsShepp.LpRmSmLm(x, -y, -phi); % reflect
            L = abs(t) + abs(u) + abs(v);
            if (valid && Lmin > L)        
                path = ReedsShepp(ReedsShepp.PATH_TYPES(6,:), t, -.5*pi, u, v);
                Lmin = L;
            end

            [valid, t, u, v] = ReedsShepp.LpRmSmLm(-x, -y, phi); % timeflip + reflect
            L = abs(t) + abs(u) + abs(v);
            if (valid && Lmin > L)        
                path = ReedsShepp(ReedsShepp.PATH_TYPES(6,:), -t, .5*pi, -u, -v);
                Lmin = L;
            end

            [valid, t, u, v] = ReedsShepp.LpRmSmRm(x, y, phi);
            L = abs(t) + abs(u) + abs(v);
            if (valid && Lmin > L)        
                path = ReedsShepp(ReedsShepp.PATH_TYPES(9,:), t, -.5*pi, u, v);
                Lmin = L;
            end

            [valid, t, u, v] = ReedsShepp.LpRmSmRm(-x, y, -phi); % timeflip
            L = abs(t) + abs(u) + abs(v);
            if (valid && Lmin > L)        
                path = ReedsShepp(ReedsShepp.PATH_TYPES(9,:), -t, .5*pi, -u, -v);
                Lmin = L;
            end

            [valid, t, u, v] = ReedsShepp.LpRmSmRm(x, -y, -phi); % reflect
            L = abs(t) + abs(u) + abs(v);
            if (valid && Lmin > L)        
                path = ReedsShepp(ReedsShepp.PATH_TYPES(10,:), t, -.5*pi, u, v);
                Lmin = L;
            end        

            [valid, t, u, v] = ReedsShepp.LpRmSmRm(-x, -y, phi); % timeflip + reflect
            L = abs(t) + abs(u) + abs(v);
            if (valid && Lmin > L)        
                path = ReedsShepp(ReedsShepp.PATH_TYPES(10,:), -t, .5*pi, -u, -v);
                Lmin = L;
            end           


            % backwards
            xb = x*cos(phi) + y*sin(phi);
            yb = x*sin(phi) - y*cos(phi);

            [valid, t, u, v] = ReedsShepp.LpRmSmLm(xb, yb, phi);
            L = abs(t) + abs(u) + abs(v);
            if (valid && Lmin > L)        
                path = ReedsShepp(ReedsShepp.PATH_TYPES(7,:), v, u, -.5*pi, t);
                Lmin = L;
            end           

            [valid, t, u, v] = ReedsShepp.LpRmSmLm(-xb, yb, -phi); % timeflip
            L = abs(t) + abs(u) + abs(v);
            if (valid && Lmin > L)        
                path = ReedsShepp(ReedsShepp.PATH_TYPES(7,:), -v, -u, .5*pi, -t);
                Lmin = L;
            end

            [valid, t, u, v] = ReedsShepp.LpRmSmLm(xb, -yb, -phi); % reflect
            L = abs(t) + abs(u) + abs(v);
            if (valid && Lmin > L)        
                path = ReedsShepp(ReedsShepp.PATH_TYPES(8,:), v, u, -.5*pi, t);
                Lmin = L;
            end

            [valid, t, u, v] = ReedsShepp.LpRmSmLm(-xb, -yb, phi); % timeflip + reflect
            L = abs(t) + abs(u) + abs(v);
            if (valid && Lmin > L)        
                path = ReedsShepp(ReedsShepp.PATH_TYPES(8,:), -v, -u, .5*pi, -t);
                Lmin = L;
            end        

            [valid, t, u, v] = ReedsShepp.LpRmSmRm(xb, yb, phi);
            L = abs(t) + abs(u) + abs(v);
            if (valid && Lmin > L)        
                path = ReedsShepp(ReedsShepp.PATH_TYPES(11,:), v, u, -.5*pi, t);
                Lmin = L;
            end        

            [valid, t, u, v] = ReedsShepp.LpRmSmRm(-xb, yb, -phi); % timeflip
            L = abs(t) + abs(u) + abs(v);
            if (valid && Lmin > L)        
                path = ReedsShepp(ReedsShepp.PATH_TYPES(11,:), -v, -u, .5*pi, -t);
                Lmin = L;
            end 

            [valid, t, u, v] = ReedsShepp.LpRmSmRm(xb, -yb, -phi); % reflect
            L = abs(t) + abs(u) + abs(v);
            if (valid && Lmin > L)        
                path = ReedsShepp(ReedsShepp.PATH_TYPES(12,:), v, u, -.5*pi, t);
                Lmin = L;
            end         

            [valid, t, u, v] = ReedsShepp.LpRmSmRm(-xb, -yb, phi); % timeflip + reflect
            L = abs(t) + abs(u) + abs(v);
            if (valid && Lmin > L)        
                path = ReedsShepp(ReedsShepp.PATH_TYPES(12,:), -v, -u, .5*pi, -t);
                Lmin = L;
            end           
        end


        % formula 8.11 *** TYPO IN PAPER ***    
        function [valid, t,u,v] = LpRmSLmRp(x, y, phi) 
            u=0; t=0; v=0;
            xi = x + sin(phi);
            eta = y - 1. - cos(phi);
            [rho, theta] = ReedsShepp.polar(xi, eta);
            if (rho >= 2.)        
                u = 4. - sqrt(rho*rho - 4.);
                if (u <= ReedsShepp.ZERO)            
                    t = ReedsShepp.mod2pi(atan2((4-u)*xi -2*eta, -2*xi + (u-4)*eta));
                    v = ReedsShepp.mod2pi(t - phi);
                    assert(abs(4*sin(t)-2*cos(t)-u*sin(t)-sin(phi) - x) < ReedsShepp.RS_EPS);
                    assert(abs(-4*cos(t)-2*sin(t)+u*cos(t)+cos(phi)+1 - y) < ReedsShepp.RS_EPS);
                    assert(abs(ReedsShepp.mod2pi(t-v-phi)) < ReedsShepp.RS_EPS);                
                    valid = (t>=-ReedsShepp.ZERO && v>=-ReedsShepp.ZERO);
                    return            
                end
            end
            valid = false;
            return;
        end

        function path = CCSCC(x, y, phi, path, backwards_enabled)    
            Lmin = path.length() - pi;

            if (~backwards_enabled)
                return; % all CCSCC paths involves backwards driving
            end            
            
            [valid, t, u, v] = ReedsShepp.LpRmSLmRp(x, y, phi);
            L = abs(t) + abs(u) + abs(v);
            if (valid && Lmin > L)        
                path = ReedsShepp(ReedsShepp.PATH_TYPES(17,:), t, -.5*pi, u, -.5*pi, v);
                Lmin = L;
            end    

            [valid, t, u, v] = ReedsShepp.LpRmSLmRp(-x, y, -phi); % timeflip
            L = abs(t) + abs(u) + abs(v);
            if (valid && Lmin > L)        
                path = ReedsShepp(ReedsShepp.PATH_TYPES(17,:), -t, .5*pi, -u, .5*pi, -v);
                Lmin = L;
            end  

            [valid, t, u, v] = ReedsShepp.LpRmSLmRp(x, -y, -phi); % reflect
            L = abs(t) + abs(u) + abs(v);
            if (valid && Lmin > L)        
                path = ReedsShepp(ReedsShepp.PATH_TYPES(18,:), t, -.5*pi, u, -.5*pi, v);
                Lmin = L;
            end          

            [valid, t, u, v] = ReedsShepp.LpRmSLmRp(-x, -y, phi); % timeflip + reflect
            L = abs(t) + abs(u) + abs(v);
            if (valid && Lmin > L)        
                path = ReedsShepp(ReedsShepp.PATH_TYPES(18,:), -t, .5*pi, -u, .5*pi, -v);
                Lmin = L;
            end           
        end
    end
end