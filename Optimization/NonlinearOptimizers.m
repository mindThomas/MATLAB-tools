classdef NonlinearOptimizers

    % Nonlinear optimization
    %
    %  min   f
    % theta
    %
    % f = nonlinear objective         [1 x 1]
    % theta = optimization variables  [n x 1]
    %
    % f, grad, and hess should be supplied as functions, e.g.:       
    % f = @(theta) ...                         [1 x 1]
    % grad = df/dtheta = @(theta) ...          [n x 1]
    % hess = (df/dtheta)/dtheta = @(theta) ... [n x n]
    methods (Static)

        function theta_opt = gradient_descent(theta0, f, grad, varargin)
            % Stopping criteria threshold
            if (isempty(varargin))
                epsilon = 0.01;
            else
                epsilon = varargin{1};
            end

            % Define the step size
            step_size = norm(theta0);
            if (step_size < eps)
                step_size = 0.01;
            end

            theta = theta0;
            while (1)
                % Steps in the negative gradient direction
                delta = -step_size * grad(theta);                

                % Stopping criteria
                if (norm(delta / theta) < epsilon)
                    break;
                end

                % Update/step in the direction
                theta = theta + delta;
            end           
            theta_opt = theta;
        end

        function theta_opt = newtons_method(theta0, f, grad, hess, varargin)
            % Approximates the neighborhood around the initial guess
            % as a quadratic function
            % f(theta0 + dt) = f(theta0) + grad(theta0) * dt + hess(theta0) * dt^2

            % Stopping criteria threshold
            if (isempty(varargin))
                epsilon = 0.01;
            else
                epsilon = varargin{1};
            end

            theta = theta0;
            while (1)               
                H = hess(theta);
                % Check eigenvalues to ensure that the approximated
                % objective function is convex
                if (min(eig(H)) >= 0)              
                    % Perform optimal step by jumping to the bottom of the
                    % qudratic approximated valley
                    delta = -H \ grad(theta); %-inv(H) * grad(theta);
                else
                    warning('Negative (semi)-definite quadratic approximation. Doing gradient descent instead.');
                    % Steps in the negative gradient direction
                    delta = -0.01 * grad(theta);
                end                

                % Stopping criteria
                if (norm(delta / theta) < epsilon)
                    break;
                end

                % Update/step in the direction
                theta = theta + delta;
            end           
            
            theta_opt = theta;
        end

    end    
    

    % Sum of squares non-linear optimization
    %
    %  min   r'*r
    % theta
    %
    % r = nonlinear residual function  [m x 1]
    % theta = optimization variables   [n x 1]
    %
    % r, J, and H should be supplied as functions, e.g.:       
    % r = @(theta) ...                         [m x 1]
    % J = jacobian = dr/dtheta = @(theta) ...  [m x n]   
    methods (Static)

        function theta_opt = gauss_newton(theta0, r, Jr, varargin)
            % Stopping criteria threshold
            if (isempty(varargin))
                epsilon = 0.01;
            else
                epsilon = varargin{1};
            end

            theta = theta0;
            while (1)
                % Perform optimal step by jumping to the bottom of the
                % qudratic valley (since objective function is quadratic)
                J = Jr(theta);
                J_pseudoinv = (J'*J) \ J'; %inv(J'*J) * J';                

                delta = -J_pseudoinv * r(theta);

                % Update/step in the direction
                theta = theta + delta;                
                
                % Stopping criteria
                if (norm(delta / (theta - theta0 + eps)) < epsilon)
                    break;
                end
            end
            
            theta_opt = theta;
        end

        function theta_opt = gauss_newton_data_fitting(theta0, f, Jf, y, varargin)
            % Gauss newton where the residual function is defined as:
            % r = y - f
            
            % Stopping criteria threshold
            if (isempty(varargin))
                epsilon = 1e-8;
            else
                epsilon = varargin{1};
            end

            theta = theta0;
            it = 0;
            while (1)
                it = it + 1;
                
                % Perform optimal step by jumping to the bottom of the
                % qudratic valley (since objective function is quadratic)
                J = Jf(theta);
                                
                if (cond(J'*J) > 1e8)
                    warning('Stopping due to non-invertible Jacobian');
                    break;
                end
                
                J_pseudoinv = (J'*J) \ J'; %inv(J'*J) * J';                

                % Compute the residual
                r = y - f(theta);
                
                delta = J_pseudoinv * r;

                % Update/step in the direction
                theta = theta + delta;                
                
                % Stopping criteria
                if (max(delta ./ (theta - theta0 + eps)) < epsilon)
                    break;
                end
            end     
            
            theta_opt = theta;
            %disp(sprintf('Iterations: %d\n', it));
        end       

    end      
   
end