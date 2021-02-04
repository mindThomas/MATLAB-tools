% B-spline with uniformly distributed knots with fixed distance of 1
% This simplifies the evaluation functions
%
% This class also includes B-spline fitting functions
%
classdef Bspline_uniform
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
        
        function kappa = curvature(obj, t)
            % Return signed curvature (positive when turning left/counter-clockwise)
            kappa = obj.curvature_(obj.p, obj.control_points, t);
        end
        
        function t = knots(obj)
            % Get the full knot vector of the current B-spline (uniformly
            % distributed)
            % Number of elements in knot vector = spline.n
            t_min = 0;
            t_max = obj.n - 1;
            t = t_min:t_max;
        end
        
        function t = effective_knots(obj)
            % Get the effective knot vector of the current B-spline (uniformly
            % distributed)
            % Number of elements in knot vector = spline.original_n
            t_min = ((obj.n-obj.original_n)/2);
            t_max = (obj.original_n-1+(obj.n-obj.original_n)/2);
            t = t_min:t_max;
        end
        
        function t = find_t(obj, point, varargin)            
            % Find progress parameter t given a query point
            % and possibly an initial guess for t
            
            knots = obj.knots();
            effective_knots = obj.effective_knots();
            t_min = min(effective_knots);
            t_max = max(effective_knots);
            
            % Define initial guess for t
            if ~isempty(varargin)
                t = varargin{1};
            else
                %t_min = (obj.n-obj.original_n)/2;
                %t_max = (obj.original_n-1+(obj.n-obj.original_n)/2);
                %t = (t_min + t_max) / 2; % midpoint as starting point
                % Find nearest control point and use its knot value as
                % initial guess
                sq_dist = diag((point - obj.control_points)' * (point - obj.control_points));
                [~,idx] = min(sq_dist);                
                t = knots(idx);
                t = min(t, t_max);
                t = max(t, t_min);
            end
            
            % Iterative algorithm to find t            
            while (1)
                % Compute curvature at current guess
                kappa = obj.curvature(t);
                R = 1/kappa; % turn radius
                
                % Compute longitudinal and lateral directions at guess
                df_dt = obj.deval_dt(obj.p, obj.control_points, t);
                lon_dir = df_dt / norm(df_dt);
                lat_dir = [0,1;-1,0] * lon_dir; % rotate so the vector points right
                
                % Compute lateral and longitudinal difference between
                % current guess location location on the B-spline and the
                % query point
                q = obj.evaluate(t);
                v = point - q;                
                lon_err = lon_dir' * v;
                lat_err = lat_dir' * v;
                
                % Compute optimal arc angle to move the guess to minimize
                % the longitudinal difference
                if (abs(lat_err) < abs(R))
                    theta = atan(lon_err / (R+lat_err));
                else
                    % Query point is more than a radius away from the guess.
                    % Just add 90 degrees in the direction that minimizes the longitudinal error
                    theta = sign(lon_err) * deg2rad(90);
                end
                
                % Compute arc curve length to move the progress using the
                % and convert that to actual delta in t based on progress stretch
                arc_curve_length_opt = R * theta;
                stretch = obj.cartesian_stretch(obj.p, obj.control_points, t);
                t = t + arc_curve_length_opt / stretch;
                t = min(t, t_max);
                t = max(t, t_min);
                
                % Termination criteria
                if (abs(theta) < deg2rad(1))
                    break;
                end        
            end            
        end
    end
    
    methods (Static)     
        
        function f = eval(p, control_points, t)
            n = length(control_points);
            f = 0;
            for i = 0:(n-1)
                f = f + control_points(:,i+1) * Bspline_uniform.basis(p, i, t);
            end
        end
        
        function df = deval_dt(p, control_points, t)
            % deval / dt
            n = length(control_points);
            df = 0;
            for i = 0:(n-1)
                df = df + control_points(:,i+1) * Bspline_uniform.dbasis(p, i, t);
            end
        end 
        
        function ddf = ddeval_dt_dt(p, control_points, t)
            % ddeval / dt / dt
            n = length(control_points);
            ddf = 0;
            for i = 0:(n-1)
                ddf = ddf + control_points(:,i+1) * Bspline_uniform.ddbasis(p, i, t);
            end
        end         

        function df = deval_dp(p, control_points, t)
            % deval / dcontrol_points
            n = length(control_points);
            df = zeros(1, n);
            for i = 0:(n-1)
                df(:, i+1) = Bspline_uniform.basis(p, i, t);
            end
        end                 
        
        function heading = heading(p, control_points, t)
            df = Bspline_uniform.deval_dt(p, control_points, t);
            heading = atan2(df(2), df(1));
        end
        
        function dheading = dheading(p, control_points, t)                        
            % d atan2 / dx = -y / (x^2 + y^2)
            % d atan2 / dy = x / (x^2 + y^2)
            % d heading / dt = dheading/df * df/dt
            df = Bspline_uniform.deval_dt(p, control_points, t);
            dheading_df = [-df(2), df(1)] / (df(1)^2 + df(2)^2);
            df_dt = Bspline_uniform.ddeval_dt_dt(p, control_points, t);
            dheading = dheading_dp * df_dt;
        end
        
        function stretch = cartesian_stretch(p, control_points, t)
            % Computes the cartesian stretch which maps between the
            % uniformly distributed knots
            % Stretch = cartesian distance / progress
            df = Bspline_uniform.deval_dt(p, control_points, t);
            stretch = norm(df);
        end       
        
        function kappa = curvature_(p, control_points, t)
            % Return signed curvature (positive when turning left/counter-clockwise)
            df = Bspline_uniform.deval_dt(p, control_points, t);
            ddf = Bspline_uniform.ddeval_dt_dt(p, control_points, t);            
            
            stretch = sqrt(df'*df); % same as cartesian_stretch
            
            kappa = (df(1)*ddf(2) - df(2)*ddf(1)) / stretch^3;
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
                B = (u - i)/p * Bspline_uniform.basis(p-1, i, s, u);
                B = B + ((i+p+1) - u) / ((i+p+1) - (i+1)) * Bspline_uniform.basis(p-1, i+1, s, u);
                %B = B + (1 - (u - (i+1)) / ((i+p+1) - (i+1))) * Bspline_uniform.basis(p-1, i+1, s, u);                     
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
                dB =      1/p * Bspline_uniform.basis(p-1, i, s, u);
                dB = dB + (u - i)/p * Bspline_uniform.dbasis(p-1, i, s, u);
                dB = dB + -1 / ((i+p+1) - (i+1)) * Bspline_uniform.basis(p-1, i+1, s, u);
                dB = dB + ((i+p+1) - u) / ((i+p+1) - (i+1)) * Bspline_uniform.dbasis(p-1, i+1, s, u);                
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
                ddB =      1/p * Bspline_uniform.dbasis(p-1, i, s, u);
                ddB = ddB + 1/p * Bspline_uniform.dbasis(p-1, i, s, u);
                ddB = ddB + (u - i)/p * Bspline_uniform.ddbasis(p-1, i, s, u);
                ddB = ddB + -1 / ((i+p+1) - (i+1)) * Bspline_uniform.dbasis(p-1, i+1, s, u);
                ddB = ddB + -1 / ((i+p+1) - (i+1)) * Bspline_uniform.dbasis(p-1, i+1, s, u);                
                ddB = ddB + ((i+p+1) - u) / ((i+p+1) - (i+1)) * Bspline_uniform.ddbasis(p-1, i+1, s, u);                
            end
        end          

        %% Reference optimization
        function obj = fit_reference(xr, yr, num_control_points, degree)            
            obj = Bspline_uniform;
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
        
        %% Objective for fitting B-spline with
        % regularizer on the distance between control points
        function cost = objective(p, xr, yr, tr, Px, Py) 
            assert(length(xr) == length(yr));
            assert(length(xr) == length(tr));
            assert(length(Px) == length(Py));            
            cost = 0;
            m = length(xr);
            n = length(Px);   
            
            mu = 0.1; % inequality weight
            
            for (j = 1:m)
                cost = cost + (xr(j) - Bspline_uniform.eval(p, Px, tr(j)))^2;
                cost = cost + (yr(j) - Bspline_uniform.eval(p, Py, tr(j)))^2;
                
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
            
            gamma = 0.5 * m/(n-2);
            normalizer = mean(sum(diff([xr', yr']).^2, 2));
            cost = cost + gamma/normalizer*sum(conv(Px, [0.5, -1, 0.5], 'valid').^2);
            cost = cost + gamma/normalizer*sum(conv(Py, [0.5, -1, 0.5], 'valid').^2);
        end
        
        function [dcost_dtr, dcost_dPx, dcost_dPy] = gradient(p, xr, yr, tr, Px, Py)            
            assert(length(xr) == length(yr));
            assert(length(xr) == length(tr));
            assert(length(Px) == length(Py));               
            m = length(xr);
            n = length(Px);            
                        
            dcost_dtr = zeros(1, length(tr));             
            dcost_dPx = zeros(1, length(Px));
            dcost_dPy = zeros(1, length(Py));            
            
            mu = 0.1; % inequality weight
            
            for (j = 1:m)
                % dcost / dtr
                dfx_dt = Bspline_uniform.deval_dt(p, Px, tr(j));
                dfy_dt = Bspline_uniform.deval_dt(p, Py, tr(j));
                                                
                dcost_dtr(j) = -2 * (xr(j) - Bspline_uniform.eval(p, Px, tr(j))) * dfx_dt ...
                               -2 * (yr(j) - Bspline_uniform.eval(p, Py, tr(j))) * dfy_dt;                

%                 if (tr(j) < 1)
%                     dcost_dtr(j) = dcost_dtr(j) - mu/(max(0 + tr(j), eps));
%                 end
%                 if (tr(j) > n-2)
%                     dcost_dtr(j) = dcost_dtr(j) + mu/(max(n-1 - tr(j), eps));
%                 end
                           
                           
                % dcost / dPx and dcost / dPy
                dfx_Px = Bspline_uniform.deval_dp(p, Px, tr(j));
                dfy_Py = Bspline_uniform.deval_dp(p, Py, tr(j));
           
                dcost_dPx = dcost_dPx + -2 * (xr(j) - Bspline_uniform.eval(p, Px, tr(j))) * dfx_Px;
                dcost_dPy = dcost_dPy + -2 * (yr(j) - Bspline_uniform.eval(p, Py, tr(j))) * dfy_Py;                       
            end                       
            

        end
        
        function dcost = gradient2(p, xr, yr, tr, Px, Py)
            [dcost_dtr, dcost_dPx, dcost_dPy] = Bspline_uniform.gradient(p, xr, yr, tr, Px, Py);
            dcost = [dcost_dtr, dcost_dPx, dcost_dPy];
        end        
        

        %% Objective for fitting B-spline with a different regularizer
        % Instead of regularizing on the distance between control points
        % we regularize the distance between the control point associated 
        % with a given knot and the corresponding queried position at the knot        
        function cost = objective_dist_regularizer(p, xr, yr, tr, Px, Py, tc) 
            assert(length(xr) == length(yr));
            assert(length(xr) == length(tr));
            assert(length(Px) == length(Py));
            %assert(length(tc) == length(Px));
            cost = 0;
            m = length(xr);
            n = length(Px);   
            
            mu = 0.1; % inequality weight
            
            for (j = 1:m)
                cost = cost + (xr(j) - Bspline_uniform.eval(p, Px, tr(j)))^2;
                cost = cost + (yr(j) - Bspline_uniform.eval(p, Py, tr(j)))^2;
                
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
                cost = cost + gamma*(Px(i) - Bspline_uniform.eval(p, Px, tc(i)))^2;
                cost = cost + gamma*(Py(i) - Bspline_uniform.eval(p, Py, tc(i)))^2;               
            end            
        end
        
        function [dcost_dtr, dcost_dPx, dcost_dPy, dcost_dtc] = gradient_dist_regularizer(p, xr, yr, tr, Px, Py, tc)            
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
                dfx_dt = Bspline_uniform.deval_dt(p, Px, tr(j));
                dfy_dt = Bspline_uniform.deval_dt(p, Py, tr(j));
                                                
                dcost_dtr(j) = -2 * (xr(j) - Bspline_uniform.eval(p, Px, tr(j))) * dfx_dt ...
                               -2 * (yr(j) - Bspline_uniform.eval(p, Py, tr(j))) * dfy_dt;                

%                 if (tr(j) < 1)
%                     dcost_dtr(j) = dcost_dtr(j) - mu/(max(0 + tr(j), eps));
%                 end
%                 if (tr(j) > n-2)
%                     dcost_dtr(j) = dcost_dtr(j) + mu/(max(n-1 - tr(j), eps));
%                 end
                           
                           
                % dcost / dPx and dcost / dPy
                dfx_Px = Bspline_uniform.deval_dp(p, Px, tr(j));
                dfy_Py = Bspline_uniform.deval_dp(p, Py, tr(j));
           
                dcost_dPx = dcost_dPx + -2 * (xr(j) - Bspline_uniform.eval(p, Px, tr(j))) * dfx_Px;
                dcost_dPy = dcost_dPy + -2 * (yr(j) - Bspline_uniform.eval(p, Py, tr(j))) * dfy_Py;                       
            end   
            
            gamma = 1;
            for (i = 1:n)
                % dcost / dtr
                dfx_dt = Bspline_uniform.deval_dt(p, Px, tc(i));
                dfy_dt = Bspline_uniform.deval_dt(p, Py, tc(i));
                                                
                dcost_dtc(i) = -2 * gamma * (Px(i) - Bspline_uniform.eval(p, Px, tc(i))) * dfx_dt ...
                               -2 * gamma * (Py(i) - Bspline_uniform.eval(p, Py, tc(i))) * dfy_dt;                   
            end            
            

        end
        
        function dcost = gradient_dist_regularizer2(p, xr, yr, tr, Px, Py, tc)
            [dcost_dtr, dcost_dPx, dcost_dPy, dcost_dtc] = Bspline_uniform.gradient_dist_regularizer(p, xr, yr, tr, Px, Py, tc);
            dcost = [dcost_dtr, dcost_dPx, dcost_dPy, dcost_dtc];
        end          

        
        %% Find locations along the B-spline based on reference points
        function cost = objective_find_location(p, xr, yr, tr, Px, Py) 
            assert(length(xr) == length(yr));
            assert(length(xr) == length(tr));
            assert(length(Px) == length(Py));               
            cost = 0;
            m = length(xr);            
            
            for (j = 1:m)
                cost = cost + (xr(j) - Bspline_uniform.eval(p, Px, tr(j)))^2;
                cost = cost + (yr(j) - Bspline_uniform.eval(p, Py, tr(j)))^2;
            end
        end
        
        function dcost_dtr = gradient_find_location(p, xr, yr, tr, Px, Py)
            assert(length(xr) == length(yr));
            assert(length(xr) == length(tr));
            assert(length(Px) == length(Py));               
            m = length(tr);                      
                        
            dcost_dtr = zeros(1, length(tr));             
            
            for (j = 1:m)
                % dcost / dtr
                dfx_dt = Bspline_uniform.deval_dt(p, Px, tr(j));
                dfy_dt = Bspline_uniform.deval_dt(p, Py, tr(j));
                                                
                dcost_dtr(j) = -2 * (xr(j) - Bspline_uniform.eval(p, Px, tr(j))) * dfx_dt ...
                               -2 * (yr(j) - Bspline_uniform.eval(p, Py, tr(j))) * dfy_dt;                                                      
            end                       
        end
        
        
        function ddcost_dtr = hessian_find_location(p, xr, yr, tr, Px, Py)
            assert(length(xr) == length(yr));
            assert(length(xr) == length(tr));
            assert(length(Px) == length(Py));               
            m = length(tr);                      
                        
            ddcost_dtr = zeros(length(tr), length(tr));             
            
            for (j = 1:m)
                % dcost / dtr
                dfx_dt = Bspline_uniform.deval_dt(p, Px, tr(j));
                dfy_dt = Bspline_uniform.deval_dt(p, Py, tr(j));
                ddfx_dt = Bspline_uniform.ddeval_dt_dt(p, Px, tr(j));
                ddfy_dt = Bspline_uniform.ddeval_dt_dt(p, Py, tr(j));                                                
                
                ddcost_dtr(j,j) =   2 * dfx_dt^2 ...
                                  - 2 * (xr(j) - dfx_dt) * ddfx_dt ...
                                  + 2 * dfy_dt^2 ...
                                  - 2 * (yr(j) - dfy_dt) * ddfy_dt;
            end                       
        end        
        
        function tr = find_locations(p, xr, yr, tr, Px, Py)
            step = 0.1;
            
            dcost_dtr = inf;
            
            while (max(abs(dcost_dtr)) > 0.01)
                dcost_dtr = Bspline_uniform.gradient_find_location(p, xr, yr, tr, Px, Py);                
                max(abs(dcost_dtr))
                tr = tr - step*dcost_dtr;                                                
            end
        end
        
        function tr = find_locations2(p, xr, yr, tr, Px, Py)                                    
            m = length(tr);    
            tr0 = tr;
            
            % Formulate problem as minimization of Sum of Squares
            % such that we can use Gauss Newton or Levenberg-Marquardt
            % f = [
            %    cost = cost + (xr(j) - Bspline_uniform.eval(p, Px, tr(j)))^2;
            %    cost = cost + (yr(j) - Bspline_uniform.eval(p, Py, tr(j)))^2;     
            %    ...
            % ]
            
            delta = inf;
            it = 0;
            while (max(delta ./ (tr + eps)) >= 0.1)
                it = it + 1;
                
                % Residual vector
                r = zeros(2*m, 1);
                for (j = 1:m)
                    r(2*(j-1)+1) = xr(j) - Bspline_uniform.eval(p, Px, tr(j));
                    r(2*(j-1)+2) = yr(j) - Bspline_uniform.eval(p, Py, tr(j));
                end

                % Jacobian of the B-spline function
                J = zeros(2*m, m);
                for (j = 1:m)
                    J(2*(j-1)+1, j) = Bspline_uniform.deval_dt(p, Px, tr(j));
                    J(2*(j-1)+2, j) = Bspline_uniform.deval_dt(p, Py, tr(j));
                end

                J_pseudoinv = (J'*J) \ J'; %inv(J'*J) * J';          
                
                delta = (J_pseudoinv * r)';
                tr = tr + delta;                
            end
            
            %disp(sprintf('Iterations: %d\n', it));
        end
        
        function tr = find_locations3(p, xr, yr, tr, Px, Py)                                    
            m = length(tr);    
            tr0 = tr;
            
            % Formulate problem as minimization of Sum of Squares
            % such that we can use Gauss Newton or Levenberg-Marquardt
            % f = [
            %    cost = cost + (xr(j) - Bspline_uniform.eval(p, Px, tr(j)))^2;
            %    cost = cost + (yr(j) - Bspline_uniform.eval(p, Py, tr(j)))^2;     
            %    ...
            % ]
            
            for (j = 1:m)                
                delta_first = inf;
                %for (i = 1:1000)
                it = 0;
                while (1)
                    it = it + 1;
                    % Residual vector
                    r = [xr(j) - Bspline_uniform.eval(p, Px, tr(j));
                         yr(j) - Bspline_uniform.eval(p, Py, tr(j))];

                    % Jacobian of the B-spline function
                    J = [Bspline_uniform.deval_dt(p, Px, tr(j));
                         Bspline_uniform.deval_dt(p, Py, tr(j))];                

                    J_pseudoinv = inv(J'*J) * J';                

                    delta = (J_pseudoinv * r)';
                    tr(j) = tr(j) + delta;
                                        
                    if (delta_first == inf)
                        delta_first = delta;
                    end
                    
%                     % Stop criteria
%                     if (delta / delta_first < 0.1)
%                         break;
%                     end
                    % Stopping criteria
                    if (max(delta ./ (tr(j) - tr0(j) + eps)) < 0.1)
                        break;
                    end
                end
                %disp(sprintf('Iterations: %d\n', it));
            end
        end   
        
        function tr = find_locations4(p, xr, yr, tr, Px, Py)
            for (i = 1:size(tr, 2))    
                tr(i) = NonlinearOptimizers.gauss_newton_data_fitting(tr(i), ...
                    @(tr) [Bspline_uniform.eval(p, Px, tr);
                           Bspline_uniform.eval(p, Py, tr)], ...
                    @(tr) [Bspline_uniform.deval_dt(p, Px, tr);
                           Bspline_uniform.deval_dt(p, Py, tr)], ...   
                          [xr(i); yr(i)], ...
                      1e-4);
            end
        end

        function tr = find_locations5(p, xr, yr, tr, Px, Py)
            f = {};
            J = {};
            for (i = 1:size(tr, 2))    
                f{end+1} = @(tr) [Bspline_uniform.eval(p, Px, tr(i));
                           Bspline_uniform.eval(p, Py, tr(i))];
                J{end+1} = @(tr) [Bspline_uniform.deval_dt(p, Px, tr(i));
                           Bspline_uniform.deval_dt(p, Py, tr(i))];
            end

            fun_f = @(x) cell2mat(cellfun(@(ff) ff(x)', f, 'UniformOutput', false))';
            J_tmp = @(x) blkdiag(cellfun(@(ff) ff(x), J, 'UniformOutput', false));
            cellblkdiag = @(cell_array) blkdiag(cell_array{:});
            fun_J = @(x)cellblkdiag(J_tmp(x));
            y_ref = [xr'; yr'];
            y_ref = y_ref(:);

            tr = NonlinearOptimizers.gauss_newton_data_fitting(tr', ...
                fun_f, ...
                fun_J, ...
                y_ref, ...
                1e-4)';
        end
        
        function tr = find_locations_fast(p, xr, yr, Px, Py)
            
        end        
        
    end
    

end