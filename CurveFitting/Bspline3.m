classdef Bspline3
    properties %(SetAccess = private)
        p % degree / order        
        n % number of control points        
        control_points   % matrix of control points stored as column vectors, rows=elements of control points, columns=different control points   
        original_n
        original_control_points
        tr % reference query locations
    end
    methods
        function f = evaluate(obj, t)
            f = obj.eval(obj.p, obj.control_points, t);
        end
    end
    methods (Static)     
        
        function f = eval(p, control_points, t)
            n = length(control_points);
            f = 0;
            for i = 0:(n-1)
                f = f + control_points(:,i+1) * Bspline3.basis(p, i, t);
            end
        end
        
        function df = deval_dt(p, control_points, t)
            % deval / dt
            n = length(control_points);
            df = 0;
            for i = 0:(n-1)
                df = df + control_points(:,i+1) * Bspline3.dbasis(p, i, t);
            end
        end 

        function df = deval_dp(p, control_points, t)
            % deval / dt
            n = length(control_points);
            df = zeros(1, n);
            for i = 0:(n-1)
                df(:, i+1) = Bspline3.basis(p, i, t);
            end
        end     
        
        
        
        
        function ddf = ddeval(p, control_points, t)
            % deval / dt
            n = length(control_points);
            ddf = 0;
            for i = 0:(n-1)
                ddf = ddf + control_points(:,i+1) * Bspline3.ddbasis(p, i, t);
            end
        end          
        
        function heading = heading(p, control_points, u)
            df = deval(p, control_points, u);
            heading = atan2(df(2), df(1));
        end
        
        function dheading = dheading(p, control_points, u)                        
            % d atan2 / dx = -y / (x^2 + y^2)
            % d atan2 / dy = x / (x^2 + y^2)
            % d heading / du = dheading/dp * dp/du
            df = deval(p, control_points, u);
            dheading_dp = [-df(2), df(1)] / (df(1)^2 + df(2)^2);
            df_du = ddeval(p, control_points, u);
            dheading = dheading_dp * df_du;

        end
        
        % Cox-de Boor recursion formula
        % Used to evaluate the basis function at a given u associated with
        % a particular control point, i, for a B-spline basis of order p
        % See:
        % - https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-basis.html
        % - https://en.wikipedia.org/wiki/B-spline#Definition
        % In this simple Bspline the knots are uniform with length 1 seperation
        function B = basis(p, i, s, varargin)   
            if (~isempty(varargin))
                u = varargin{1};
            else
                u = s + (p+1)/2;
            end
            
            if (p == 0)                                 
                B = double(u >= i && u < i+1);                
            else
                B = (u - i)/p * Bspline3.basis(p-1, i, s, u);
                B = B + ((i+p+1) - u) / ((i+p+1) - (i+1)) * Bspline3.basis(p-1, i+1, s, u);
                %B = B + (1 - (u - (i+1)) / ((i+p+1) - (i+1))) * Bspline3.basis(p-1, i+1, s, u);                     
                % These two lines gives the same, but the latter
                % matches the definition from Wikipedia                
            end
        end             
        
        function dB = dbasis(p, i, s, varargin)      
            % dBasis / du
            if (~isempty(varargin))
                u = varargin{1};
            else
                u = s + (p+1)/2;
            end
            
            if (p == 0)                                 
                dB = 0;                
            else
                dB =      1/p * Bspline3.basis(p-1, i, s, u);
                dB = dB + (u - i)/p * Bspline3.dbasis(p-1, i, s, u);
                dB = dB + -1 / ((i+p+1) - (i+1)) * Bspline3.basis(p-1, i+1, s, u);
                dB = dB + ((i+p+1) - u) / ((i+p+1) - (i+1)) * Bspline3.dbasis(p-1, i+1, s, u);                
            end
        end             
        
        function ddB = ddbasis(p, i, s, varargin)          
            % ddBasis / du
            if (~isempty(varargin))
                u = varargin{1};
            else
                u = s + (p+1)/2;
            end
            
            if (p == 0)                                 
                ddB = 0;                
            else
                ddB =      1/p * Bspline3.dbasis(p-1, i, s, u);
                ddB = ddB + 1/p * Bspline3.dbasis(p-1, i, s, u);
                ddB = ddB + (u - i)/p * Bspline3.ddbasis(p-1, i, s, u);
                ddB = ddB + -1 / ((i+p+1) - (i+1)) * Bspline3.dbasis(p-1, i+1, s, u);
                ddB = ddB + -1 / ((i+p+1) - (i+1)) * Bspline3.dbasis(p-1, i+1, s, u);                
                ddB = ddB + ((i+p+1) - u) / ((i+p+1) - (i+1)) * Bspline3.ddbasis(p-1, i+1, s, u);                
            end
        end          

        %% Reference optimization
        function obj = fit_reference(xr, yr, num_control_points, degree)            
            obj = Bspline3;
            n = num_control_points;
            m = length(xr);            
            p = degree; % degree               
            % knots => t_i = i
            
            % Initialize reference query locations uniformly and
            j = 0:(m-1);
            tr = j*(n-1)/(m-1);
            
            % Initialize control points from equally distant reference
            % points, assuming that the reference points are approximately
            % uniformly distributed/sampled
            control_points = zeros(2, num_control_points);            
            for (i = 1:num_control_points)
                j = round((i-1)*(m-1)/(n-1)) + 1;
                control_points(:,i) = [xr(j); yr(j)];  
                
                % Update initial reference locations for those reference
                % points used to construct initial control point locations
                tr(j) = i-1;
            end

            obj.p = p;
            obj.original_n = n;
            obj.original_control_points = control_points;
            obj.n = n;
            obj.control_points = control_points;
            obj.tr = tr;
            
            % Extend the control points with an extra point at the same
            % location in the beginning and end
            obj.n = n + 2;
            obj.control_points = [control_points(:,1), control_points, control_points(:,end)]; 
            obj.tr = tr + 1;
        end
        
        function cost = objective(p, xr, yr, tr, Px, Py, tc) 
            assert(length(xr) == length(yr));
            assert(length(xr) == length(tr));
            assert(length(Px) == length(Py));
            %assert(length(tc) == length(Px));
            cost = 0;
            m = length(xr);
            n = length(Px);   
            
            mu = 0.1; % inequality weight
            
            for (j = 1:m)
                cost = cost + (xr(j) - Bspline3.eval(p, Px, tr(j)))^2;
                cost = cost + (yr(j) - Bspline3.eval(p, Py, tr(j)))^2;
                
                % Add inequality constraints
                % x < b  -->                -log(b - x)
                % x > b  -->  -x < -b  -->  -log(-b + x)
                %
%                 if (tr(j) < 1)
%                     % tr_j >= 0  -->  tr_j > 0
%                     cost = cost - mu*log(max(0 + tr(j), eps));
%                 end
%                 if (tr(j) > n-2)
%                     % tr_j < n-1
%                     cost = cost - mu*log(max(n-1 - tr(j), eps));
%                 end
            end
            
            gamma = 1;
            for (i = 1:n)
                cost = cost + gamma*(Px(i) - Bspline3.eval(p, Px, tc(i)))^2;
                cost = cost + gamma*(Py(i) - Bspline3.eval(p, Py, tc(i)))^2;               
            end            
        end
        
        function [dcost_dtr, dcost_dPx, dcost_dPy, dcost_dtc] = gradient(p, xr, yr, tr, Px, Py, tc)            
            assert(length(xr) == length(yr));
            assert(length(xr) == length(tr));
            assert(length(Px) == length(Py));   
            %assert(length(tc) == length(Px));
            m = length(xr);
            n = length(Px);            
                        
            dcost_dtr = zeros(1, length(tr));             
            dcost_dPx = zeros(1, length(Px));
            dcost_dPy = zeros(1, length(Py));
            dcost_dtc = zeros(1, length(tc)); 
            
            mu = 0.1; % inequality weight
            
            for (j = 1:m)
                % dcost / dtr
                dfx_dt = Bspline3.deval_dt(p, Px, tr(j));
                dfy_dt = Bspline3.deval_dt(p, Py, tr(j));
                                                
                dcost_dtr(j) = -2 * (xr(j) - Bspline3.eval(p, Px, tr(j))) * dfx_dt ...
                               -2 * (yr(j) - Bspline3.eval(p, Py, tr(j))) * dfy_dt;                

%                 if (tr(j) < 1)
%                     dcost_dtr(j) = dcost_dtr(j) - mu/(max(0 + tr(j), eps));
%                 end
%                 if (tr(j) > n-2)
%                     dcost_dtr(j) = dcost_dtr(j) + mu/(max(n-1 - tr(j), eps));
%                 end
                           
                           
                % dcost / dPx and dcost / dPy
                dfx_Px = Bspline3.deval_dp(p, Px, tr(j));
                dfy_Py = Bspline3.deval_dp(p, Py, tr(j));
           
                dcost_dPx = dcost_dPx + -2 * (xr(j) - Bspline3.eval(p, Px, tr(j))) * dfx_Px;
                dcost_dPy = dcost_dPy + -2 * (yr(j) - Bspline3.eval(p, Py, tr(j))) * dfy_Py;                       
            end   
            
            gamma = 1;
            for (i = 1:n)
                % dcost / dtr
                dfx_dt = Bspline3.deval_dt(p, Px, tc(i));
                dfy_dt = Bspline3.deval_dt(p, Py, tc(i));
                                                
                dcost_dtc(i) = -2 * gamma * (Px(i) - Bspline3.eval(p, Px, tc(i))) * dfx_dt ...
                               -2 * gamma * (Py(i) - Bspline3.eval(p, Py, tc(i))) * dfy_dt;                   
            end            
            

        end
        
        function dcost = gradient2(p, xr, yr, tr, Px, Py, tc)
            [dcost_dtr, dcost_dPx, dcost_dPy, dcost_dtc] = Bspline3.gradient(p, xr, yr, tr, Px, Py, tc);
            dcost = [dcost_dtr, dcost_dPx, dcost_dPy, dcost_dtc];
        end        
        
    end
    

end