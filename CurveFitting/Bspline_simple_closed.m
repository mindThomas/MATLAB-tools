classdef Bspline_simple_closed
    properties %(SetAccess = private)
        n % number of control points
        p % degree / order
        control_points   % matrix of control points stored as column vectors, rows=elements of control points, columns=different control points   
    end
    methods
        function obj = Bspline_simple_closed(control_points, degree)            
            obj.control_points = control_points; 
            obj.n = length(obj.control_points);            
            obj.p = degree; % degree            
        end        
        
        function y = eval(obj, u)
            y = 0;
            for i = 0:(obj.n-1)
                y = y + obj.control_points(:,i+1) * obj.basis(obj.p, i, mod(u+(obj.p+1)/2, obj.n)); % add the degree to center the kernel
            end
        end
        
        function dy = deval(obj, u)
            % deval / du
            dy = 0;
            for i = 0:(obj.n-1)
                dy = dy + obj.control_points(:,i+1) * obj.dbasis(obj.p, i, mod(u+(obj.p+1)/2, obj.n)); % add the degree to center the kernel                
            end
        end     
        
        function ddy = ddeval(obj, u)
            % deval / du
            ddy = 0;
            for i = 0:(obj.n-1)
                ddy = ddy + obj.control_points(:,i+1) * obj.ddbasis(obj.p, i, mod(u+(obj.p+1)/2, obj.n)); % add the degree to center the kernel                
            end
        end          
        
        function heading = heading(obj, u)
            dy = obj.deval(u);
            heading = atan2(dy(2), dy(1));
        end
        
        function dheading = dheading(obj, u)                        
            % d atan2 / dx = -y / (x^2 + y^2)
            % d atan2 / dy = x / (x^2 + y^2)
            % d heading / du = dheading/dp * dp/du
            dp = obj.deval(u);
            dheading_dp = [-dp(2), dp(1)] / (dp(1)^2 + dp(2)^2);
            dp_du = obj.ddeval(u);
            dheading = dheading_dp * dp_du;

        end
        
        % Cox-de Boor recursion formula
        % Used to evaluate the basis function at a given u associated with
        % a particular control point, i, for a B-spline basis of order p
        % See:
        % - https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-basis.html
        % - https://en.wikipedia.org/wiki/B-spline#Definition
        % In this simple Bspline the knots are uniform with length 1 seperation
        function B = basis(obj, p, i, u)            
            if (p == 0)                                 
                B = double(mod(u - i, obj.n) < 1);                
            else
                B = mod(u - i, obj.n) / mod(mod(i+p, obj.n) - i, obj.n) * obj.basis(p-1, i, u);

                %B = B + (obj.knots(i+p+2) - u) / (obj.knots(i+p+2) - obj.knots(i+2)) * obj.basis(p-1, i+1, u);
                B = B + (1 - mod(u - (i+1), obj.n) / mod(mod(i+p+1, obj.n) - mod(i+1, obj.n), obj.n)) * obj.basis(p-1, mod(i+1, obj.n), u);                     
                % These two lines gives the same, but the latter
                % matches the definition from Wikipedia                
            end
        end             
        
        function dB = dbasis(obj, p, i, u)            
            % dBasis / du
            if (p == 0)                                 
                dB = 0;
            else
                dB =      1 / mod(mod(i+p, obj.n) - i, obj.n) * obj.basis(p-1, i, u);
                dB = dB + mod(u - i, obj.n) / mod(mod(i+p, obj.n) - i, obj.n) * obj.dbasis(p-1, i, u);

                %B = B + (obj.knots(i+p+2) - u) / (obj.knots(i+p+2) - obj.knots(i+2)) * obj.basis(p-1, i+1, u);
                dB = dB + -1 / mod(mod(i+p+1, obj.n) - mod(i+1, obj.n), obj.n) * obj.basis(p-1, mod(i+1, obj.n), u);                     
                dB = dB + (1 - mod(u - (i+1), obj.n) / mod(mod(i+p+1, obj.n) - mod(i+1, obj.n), obj.n)) * obj.dbasis(p-1, mod(i+1, obj.n), u);                     
                % These two lines gives the same, but the latter
                % matches the definition from Wikipedia                
            end
        end             
        
        function ddB = ddbasis(obj, p, i, u)            
            % ddBasis / du
            if (p == 0)                                 
                ddB = 0;
            else
                ddB =       1 / mod(mod(i+p, obj.n) - i, obj.n) * obj.dbasis(p-1, i, u);
                ddB = ddB + 1 / mod(mod(i+p, obj.n) - i, obj.n) * obj.dbasis(p-1, i, u);
                ddB = ddB + mod(u - i, obj.n) / mod(mod(i+p, obj.n) - i, obj.n) * obj.ddbasis(p-1, i, u);

                %B = B + (obj.knots(i+p+2) - u) / (obj.knots(i+p+2) - obj.knots(i+2)) * obj.basis(p-1, i+1, u);
                ddB = ddB + -1 / mod(mod(i+p+1, obj.n) - mod(i+1, obj.n), obj.n) * obj.dbasis(p-1, mod(i+1, obj.n), u);                                                         
                ddB = ddB + -1 / mod(mod(i+p+1, obj.n) - mod(i+1, obj.n), obj.n) * obj.dbasis(p-1, mod(i+1, obj.n), u);                     
                ddB = ddB + (1 - mod(u - (i+1), obj.n) / mod(mod(i+p+1, obj.n) - mod(i+1, obj.n), obj.n)) * obj.ddbasis(p-1, mod(i+1, obj.n), u);                     
                % These two lines gives the same, but the latter
                % matches the definition from Wikipedia                
            end
        end          
        
    end      

end