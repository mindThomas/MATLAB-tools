% x and y should be column vectors; m x 1
% modelfun should be on the form: modelfun = @(coeffs,x) x(:, 1) ...
function coeffs = fitNonlinearModel(x, y, modelfun, init_coeffs, varargin)
    tbl = table(x, y);
    if (length(varargin) == 1) % Weighted Nonlinear Regression
        weights = varargin{1};
        mdl = fitnlm(tbl, modelfun, init_coeffs, 'Weight', weights);        
    else % Nonlinear Regression
        mdl = fitnlm(tbl, modelfun, init_coeffs);
    end
    coeffs = mdl.Coefficients{:, 'Estimate'};