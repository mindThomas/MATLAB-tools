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
        
    if (mean(diff(y)) < 0) % exponential decay
        b = y(1);
    else % assumed positive exponential step
        b = y(end);
    end
    
    % Compute initialization point of c from the measurements, as long as x
    % is increasing
    %y1 / y2 = b*exp(c*x1) / b*exp(c*x2)
    %y1 / y2 = exp(c*x1) / exp(c*x2)
    %y1 / y2 = exp(c*(x1-x2))
    %ln(y1 / y2) = c * (x1-x2)
    %c = ln(y1 / y2) / (x1-x2)
    
    x_sampled = x(1:round(length(y)/20):end);
    y_sampled = y(1:round(length(y)/20):end);    
    
    dy = y_sampled(end) - y_sampled(1);    
    idx_include = find(abs(diff(y_sampled)) > 0.1 * abs(dy));
    idx_include = unique([idx_include; idx_include+1]);
    if (mod(length(idx_include), 2) == 1)
        idx_include = idx_include(1:end-1);
    end
    if (length(idx_include) < 2)
        return; % error, not sufficient samples to estimate initialization points
    end
    x_sampled = x_sampled(idx_include);
    y_sampled = y_sampled(idx_include);
    
    c_init = mean(log(y_sampled(1:end-1) ./ y_sampled(2:end)) ./ (x_sampled(1:end-1)-x_sampled(2:end)));       
    if (isinf(c_init))
        return; % error - two x consecutive x-values have the same value
    end
    c_init = -abs(c_init);
    c = c_init;

    coeffs = [a, b, c]';
    
    delta = single([inf, inf, inf])';
    prev_delta = single([0, 0, 0])';
    new_cost = cost(coeffs(1), coeffs(2), coeffs(3));    
    while (true)
        Jr_ = Jr(coeffs(1), coeffs(2), coeffs(3));        
        delta = inv(Jr_'*Jr_) * (Jr_' * error(coeffs(1), coeffs(2), coeffs(3)));
        coeffs = coeffs - delta
        if (coeffs(3) > 0) % force negative c since we are fitting an exponentially stable function
            coeffs(3) = c_init;
        end
        
        % check for convergence
        if (norm(delta) < 1/100*abs(coeffs(2)) )
            break;
        end
        
        new_cost = cost(coeffs(1), coeffs(2), coeffs(3));
        prev_delta = delta;
    end   
    
    a = coeffs(1);
    c = coeffs(3);
    b = coeffs(2) * exp(c * (-x0));     % correct parameter with offset shift