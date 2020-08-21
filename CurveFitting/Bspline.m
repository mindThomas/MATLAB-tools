classdef Bspline 
    properties %(SetAccess = private)
        n % number of knots
        p % degree / order
        knots            % vector of knot locations (length should be the same as number of columns in control_points)
        control_points   % matrix of control points stored as column vectors, rows=elements of control points, columns=different control points   
    end
    methods
        function obj = Bspline(control_points, degree)            
            obj.control_points = control_points; 
            obj.n = length(obj.control_points);
            obj.knots = 0:(obj.n-1); % make uniform B-spline currently
            obj.p = degree; % degree
            
            % Note that since a single B-spline (basis function) already
            % spans over 1+p knots, it follows that the internal knots need to be extended with p-1 end points on each side,
            % to give full support to the first and last B-spline which affect the internal knot intervals. 
            if (obj.p > 1)
                for (i = 1:obj.p)
                    obj.knots = [obj.knots(1) - (obj.knots(2)-obj.knots(1)), obj.knots, obj.knots(end) + (obj.knots(end)-obj.knots(end-1))];
                    obj.control_points = [obj.control_points(:,1), obj.control_points, obj.control_points(:,end)];
                end
            end
            obj.n = length(obj.control_points);
        end        
        
        function y = eval(obj, u)
            y = 0;
            for i = 0:(obj.n-1)
                y = y + obj.control_points(:,i+1) * obj.basis(obj.p, i, u+(obj.p+1)/2); % add the degree to center the kernel
            end
        end
        
        % Cox-de Boor recursion formula
        % Used to evaluate the basis function at a given u associated with
        % a particular control point, i, for a B-spline basis of order p
        % See:
        % - https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-basis.html
        % - https://en.wikipedia.org/wiki/B-spline#Definition
        function B = basis(obj, p, i, u)            
            B = 0;
            if (i < 0 || i >= obj.n)
                return;
            end            
            
            if (p == 0)
                if (i+1 < obj.n)                    
                    B = double(u >= obj.knots(i+1) && u < obj.knots(i+2));
                end
            else
                if (i+p < obj.n)
                    B = B + (u - obj.knots(i+1)) / (obj.knots(i+p+1) - obj.knots(i+1)) * obj.basis(p-1, i, u);
                end
                if (i+p+1 < obj.n)
                    %B = B + (obj.knots(i+p+2) - u) / (obj.knots(i+p+2) - obj.knots(i+2)) * obj.basis(p-1, i+1, u);
                    B = B + (1 - (u - obj.knots(i+2)) / (obj.knots(i+p+2) - obj.knots(i+2))) * obj.basis(p-1, i+1, u);                     
                    % These two lines gives the same, but the latter
                    % matches the definition from Wikipedia
                end                                
            end
        end                
        
        % Closed form expressions for the basis for certain orders
        function B = basis0(obj, i, u)       
            % Check that input is valid
            % There should be 2 valid knots to query
            if (i < 0 || i+1 >= obj.n)                    
                % Error, outside of bounds
                B = 0;
                return;
            end
            
            % B = 1   if  'u' is between current and next knot
            % B = 0   otherwise
            B = double(u >= obj.knots(i+1) && u < obj.knots(i+2));            
        end
        
        function B = basis1(obj, i, u)       
            % Check that input is valid
            % There should be 3 valid knots to query
            if (i < 0 || i+2 >= obj.n)                    
                % Error, outside of bounds
                B = 0;
                return;
            end
            
            % B = ramp up     if  'u' is between knot i and knot i+1
            % B = ramp down   if  'u' is between knot i+1 and knot i+2
            % B = 0   otherwise
            if (u >= obj.knots(i+1) && u < obj.knots(i+2))
                B = (u - obj.knots(i+1)) / (obj.knots(i+2) - obj.knots(i+1));
            elseif (u >= obj.knots(i+2) && u < obj.knots(i+3))
                B = -(u - obj.knots(i+3)) / (obj.knots(i+3) - obj.knots(i+2));
            else
                B = 0;
            end
        end    
        
        function B = basis2(obj, i, u)       
            % Check that input is valid
            % There should be 4 valid knots to query
            if (i < 0 || i+3 >= obj.n)                    
                % Error, outside of bounds
                B = 0;
                return;
            end
                       
            % Based on Symbolic_CoxDeBoorRecursion.m
            %  (((s12*(u - u1))/(u1 - u2) - (s23*(u - u3))/(u2 - u3))*(u - u1))/(u1 - u3) 
            % -(((s23*(u - u2))/(u2 - u3) - (s34*(u - u4))/(u3 - u4))*(u - u4))/(u2 - u4)            
            if (u >= obj.knots(i+1) && u < obj.knots(i+2))                
                %B = (u - u1))/(u1 - u2)  *  (u - u1))/(u1 - u3);
                B = (u - obj.knots(i+1))/(obj.knots(i+1) - obj.knots(i+2))  *  (u - obj.knots(i+1))/(obj.knots(i+1) - obj.knots(i+3));
            elseif (u >= obj.knots(i+2) && u < obj.knots(i+3))
                %B = -(u - u3)/(u2 - u3)  *  (u - u1)/(u1 - u3) ...
                %    -(u - u2)/(u2 - u3)  *  (u - u4)/(u2 - u4);
                B = -(u - obj.knots(i+3))/(obj.knots(i+2) - obj.knots(i+3))  *  (u - obj.knots(i+1))/(obj.knots(i+1) - obj.knots(i+3)) ...
                    -(u - obj.knots(i+2))/(obj.knots(i+2) - obj.knots(i+3))  *  (u - obj.knots(i+4))/(obj.knots(i+2) - obj.knots(i+4));
            elseif (u >= obj.knots(i+3) && u < obj.knots(i+4))
                %B = (u - u4)/(u3 - u4)  *  (u - u4)/(u2 - u4)
                B = (u - obj.knots(i+4))/(obj.knots(i+3) - obj.knots(i+4))  *  (u - obj.knots(i+4))/(obj.knots(i+2) - obj.knots(i+4));
            else
                B = 0;
            end
        end         

        function B = basis3(obj, i, u)       
            % Check that input is valid
            % There should be 5 valid knots to query
            if (i < 0 || i+4 >= obj.n)                    
                % Error, outside of bounds
                B = 0;
                return;
            end
                       
            % Based on Symbolic_CoxDeBoorRecursion.m
            %((...
            %    (((s23*(u - u2))/(u2 - u3) - (s34*(u - u4))/(u3 - u4))*(u - u2))/(u2 - u4)...
            %   -(((s34*(u - u3))/(u3 - u4) - (s45*(u - u5))/(u4 - u5))*(u - u5))/(u3 - u5)...
            %)*(u - u5))/(u2 - u5) ...
            %-((...
            %    (((s12*(u - u1))/(u1 - u2) - (s23*(u - u3))/(u2 - u3))*(u - u1))/(u1 - u3)...
            %   -(((s23*(u - u2))/(u2 - u3) - (s34*(u - u4))/(u3 - u4))*(u - u4))/(u2 - u4)
            %)*(u - u1))/(u1 - u4)          
            if (u >= obj.knots(i+1) && u < obj.knots(i+2))     % s12              
                %B = -(u - u1)/(u1 - u2) * (u - u1)/(u1 - u3) * (u - u1)/(u1 - u4);
                B = -(u - obj.knots(i+1))/(obj.knots(i+1) - obj.knots(i+2)) * (u - obj.knots(i+1))/(obj.knots(i+1) - obj.knots(i+3)) * (u - obj.knots(i+1))/(obj.knots(i+1) - obj.knots(i+4));
            elseif (u >= obj.knots(i+2) && u < obj.knots(i+3)) % s23
                %B  = (u - u2)/(u2 - u3) * (u - u2)/(u2 - u4) * (u - u5)/(u2 - u5);
                %B += (u - u3)/(u2 - u3) * (u - u1)/(u1 - u3) * (u - u1)/(u1 - u4);
                %B += (u - u2)/(u2 - u3) * (u - u4)/(u2 - u4) * (u - u1)/(u1 - u4);
                B  = (u - obj.knots(i+2))/(obj.knots(i+2) - obj.knots(i+3)) * (u - obj.knots(i+2))/(obj.knots(i+2) - obj.knots(i+4)) * (u - obj.knots(i+5))/(obj.knots(i+2) - obj.knots(i+5)) ...
                   + (u - obj.knots(i+3))/(obj.knots(i+2) - obj.knots(i+3)) * (u - obj.knots(i+1))/(obj.knots(i+1) - obj.knots(i+3)) * (u - obj.knots(i+1))/(obj.knots(i+1) - obj.knots(i+4)) ...
                   + (u - obj.knots(i+2))/(obj.knots(i+2) - obj.knots(i+3)) * (u - obj.knots(i+4))/(obj.knots(i+2) - obj.knots(i+4)) * (u - obj.knots(i+1))/(obj.knots(i+1) - obj.knots(i+4));                
            elseif (u >= obj.knots(i+3) && u < obj.knots(i+4)) % s34
                %B  = -(u - u4)/(u3 - u4) * (u - u2)/(u2 - u4) * (u - u5)/(u2 - u5);
                %B += -(u - u3)/(u3 - u4) * (u - u5)/(u3 - u5) * (u - u5)/(u2 - u5);
                %B += -(u - u4)/(u3 - u4) * (u - u4)/(u2 - u4) * (u - u1)/(u1 - u4);
                B  = -(u - obj.knots(i+4))/(obj.knots(i+3) - obj.knots(i+4)) * (u - obj.knots(i+2))/(obj.knots(i+2) - obj.knots(i+4)) * (u - obj.knots(i+5))/(obj.knots(i+2) - obj.knots(i+5)) ...
                   + -(u - obj.knots(i+3))/(obj.knots(i+3) - obj.knots(i+4)) * (u - obj.knots(i+5))/(obj.knots(i+3) - obj.knots(i+5)) * (u - obj.knots(i+5))/(obj.knots(i+2) - obj.knots(i+5)) ...
                   + -(u - obj.knots(i+4))/(obj.knots(i+3) - obj.knots(i+4)) * (u - obj.knots(i+4))/(obj.knots(i+2) - obj.knots(i+4)) * (u - obj.knots(i+1))/(obj.knots(i+1) - obj.knots(i+4));
            elseif (u >= obj.knots(i+4) && u < obj.knots(i+5)) % s45
                %B = (u - u5)/(u4 - u5) * (u - u5)/(u3 - u5) * (u - u5)/(u2 - u5);
                B = (u - obj.knots(i+5))/(obj.knots(i+4) - obj.knots(i+5)) * (u - obj.knots(i+5))/(obj.knots(i+3) - obj.knots(i+5)) * (u - obj.knots(i+5))/(obj.knots(i+2) - obj.knots(i+5));
            else
                B = 0;
            end
        end         
    end
    

end