% x and y should be column vectors; m x 1
% modelfun should be on the form: modelfun = @(coeffs,x) x(:, 1) ...
function coeffs = fitNonlinearModel(x, y, modelfun, init_coeffs)
    tbl = table(x, y);
    mdl = fitnlm(tbl, modelfun, init_coeffs);
    coeffs = mdl.Coefficients{:, 'Estimate'};