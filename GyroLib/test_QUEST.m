function test_QUEST()

%% Create two random vectors in "body" frame
fb = randn(3,1);
mb = randn(3,1);
fb = fb./norm(fb);
mb = mb./norm(mb);
fprintf('Measurements in b-frame:\n');
fprintf('ab=[%f, %f %f]\n', fb(1), fb(2), fb(3));
fprintf('mb=[%f, %f %f]\n', mb(1), mb(2), mb(3));

%% Create random DCM from the "body" to the "navigation" frame
rz = randn;  %Z axis
ry = randn;  %Y axis
rx = randn;  %X axis
Cbn = angle_dcm(rz, ry, rx);
fprintf('\nReference DCM:\n');
disp(Cbn);

%% Rotate vectors to the "navigation" frame
fn = Cbn*fb;
mn = Cbn*mb;
fprintf('Vectors in n-frame:\n');
fprintf('an=[%f, %f %f]\n', fn(1), fn(2), fn(3));
fprintf('mn=[%f, %f %f]\n', mn(1), mn(2), mn(3));

%% Sensors weights
wf = 1.0;
wm = 1.0;

%% Call QUEST
q = QUEST(fb, mb, fn, mn, wf, wm);
Cbn_hat = quat_dcm(q);
fprintf('\nEstimated DCM:\n');
disp(Cbn_hat);

%% Check the results 
fprintf('Test - should give the identity matrix:\n');
dCnb = Cbn'*Cbn_hat;
disp(dCnb);

%% Estimation errors
[dz, dy, dx] = dcm_angle(dCnb);
fprintf('Errors = [%10.10f, %10.10f %10.10f] deg.\n', rad2deg(dx), rad2deg(dy), rad2deg(dz));

end