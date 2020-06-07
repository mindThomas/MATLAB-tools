% x and y should be column vectors; m x 1
% f(x) = a + b * exp(-c*x)
function [a, b, c] = fitExponentialDecay(x, y, varargin)
    coeffs0 = [0, y(1), 1]; % initialization point
    model = @(coeff,t)(coeff(1) + coeff(2) * exp(-coeff(3)*t(:, 1)));
    
    if (length(varargin) > 0)
        weights = varargin{1};        
        coeffs = fitNonlinearModel(x, y, model, coeffs0, 'weights', weights);
    else
        coeffs = fitNonlinearModel(x, y, model, coeffs0);
    end
    
    a = coeffs(1);
    b = coeffs(2);
    c = coeffs(3);   