x0 = [0; 0];

%% Without Jacobian
f_min = @(params)  sqrt((params(1) - 100)^2 + (params(2) - 50)^2)

opt.Jacobian='romberg';
opt.Display='iter';
[x1,~,fval_1,~,extra_arguments] = LevenbergMarquardt(f_min,x0,[],[],opt);
x1

%% With Jacobian
f_min = @(params)  sqrt((params(1) - 100)^2 + (params(2) - 50)^2)
F_min = @(params)  [ (params(1)-100) / f_min(params) ,...
                     (params(2)-50) / f_min(params) ]                        
f_combined = @(params) deal(f_min(params), F_min(params))  % combine into function with multiple outputs

opt.Jacobian='on';
opt.Display='iter';
[x1,~,fval_1,~,extra_arguments] = LevenbergMarquardt(f_combined,x0,[],[],opt);
x1