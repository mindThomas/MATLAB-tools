% x and y should be column vectors; m x 1
% f(x) = a + b * exp(c*x)
function [a, b, c] = fitExponentialDecay2(x, y) %#codegen
   
    % To make initialization point valid, subtract starting point from x
    % We can always shift the results
    x0 = x(1);
    x = x - x0;   
    n = length(x);
    
    % Formulate error function
    model = @(x, a, b, c) a + b.*exp(c.*x);
    error = @(a, b, c) y - model(x, a, b, c);
    de_dC = @(a, b, c) -[single(ones(n,1)), exp(c*x), b .* x.*exp(c*x)];    
    cost = @(a, b, c) sum( error(a, b, c).^2 );
    
    % Compute Jacobian   
    J = @(a, b, c) 2*sum( error(a, b, c) .* de_dC(a, b, c) );
    Jr = @(a, b, c) de_dC(a, b, c); % Jacobian of residual
    
    %% Perform Gauss-Newton algorithm    
    % Initialization point
    a = single(0);
    b = y(1);
    c = single(-0.1);
    
    coeffs = [a, b, c]';
    
    delta = single([inf, inf, inf])';
    prev_delta = single([0, 0, 0])';
    new_cost = cost(coeffs(1), coeffs(2), coeffs(3));    
    while (true)
        Jr_ = Jr(coeffs(1), coeffs(2), coeffs(3));        
        delta = inv(Jr_'*Jr_) * (Jr_' * error(coeffs(1), coeffs(2), coeffs(3)));
        coeffs = coeffs - delta;
        if (coeffs(3) > 0)
            coeffs(3) = single(-0.1);
        end
        
        % check for convergence
        if (norm(delta) < 1/100*abs(coeffs(3)) )
            break;
        end
        
        new_cost = cost(coeffs(1), coeffs(2), coeffs(3));
        prev_delta = delta;
    end   
    
    a = coeffs(1);
    c = coeffs(3);
    b = coeffs(2) * exp(c * (-x0));     % correct parameter with offset shift