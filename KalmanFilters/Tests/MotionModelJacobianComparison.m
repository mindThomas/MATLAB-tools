scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, '../'));
addpath(fullfile(scriptDir, '../MotionModels'));

ts = 0.05;
[f, Fx, Fu, Fq] = CoordinatedTurnModel_Discrete(ts);

% x = [ x, y, v, phi, omega ]
x = [0, 2, 0.5, deg2rad(50), deg2rad(10)]';
u = [];
q = [0, 0]';

[dfdx, dfdu, dfdq] = numjacobian3(f, x, u, q);

Fx(x, u, q)
dfdx