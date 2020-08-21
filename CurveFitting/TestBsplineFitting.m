%% Test basis function
%u = min(spline.knots):0.02:max(spline.knots);
u = 0:0.02:10;
B = zeros(size(u));
dB = zeros(size(u));
ddB = zeros(size(u));
for (i = 1:length(u))    
    B(i) = Bspline2.basis(3, 5, u(i));    
    dB(i) = Bspline2.dbasis(3, 5, u(i));    
    ddB(i) = Bspline2.ddbasis(3, 5, u(i));       
end

figure(1);
plot(u, B);
hold on;
plot(u, dB);
plot(u, ddB);
hold off;

legend('3rd order basis', '1st derivative', '2nd derivative');

%%
num_ref_points = 10;
num_control_points = 4;
order = 3;

figure(1);
xlim([0, 10]);
ylim([0, 10]);
ref = ginput(num_ref_points);
plot(ref(:,1), ref(:,2), 'x');
xlim([0, 10]);
ylim([0, 10]);

%% Initial for fitting
spline = Bspline2.fit_reference(ref(:,1), ref(:,2), num_control_points, order);

t = ((spline.n-spline.original_n)/2):0.02:(spline.original_n-1+(spline.n-spline.original_n)/2);
%t = 0:0.02:(spline.n-1);

pos = zeros(2, length(t));
for (i = 1:length(t))    
    pos(:,i) = spline.evaluate(t(i));
end

evaluation_points =  zeros(2, length(spline.tr));
for (i = 1:length(spline.tr))    
    evaluation_points(:,i) = spline.evaluate(spline.tr(i));
end

figure(1);
plot(spline.control_points(1,:), spline.control_points(2,:), 'o');
hold on;
plot(ref(:,1), ref(:,2), 'x');
plot(evaluation_points(1,:), evaluation_points(2,:), 'g^');
for (i = 1:length(spline.tr))   
    line([ref(i,1), evaluation_points(1,i)], [ref(i,2), evaluation_points(2,i)], 'Color', 'g', 'LineStyle', '--');
end
plot(pos(1,:), pos(2,:));
hold off;
axis equal;
xlim([2, 8]);
ylim([0, 10]);

%%
first = true;
while (true)    
cost = Bspline2.objective(order, ref(:,1), ref(:,2), spline.tr, spline.control_points(1,:), spline.control_points(2,:))
[dcost_dtr, dcost_dPx, dcost_dPy] = Bspline2.gradient(order, ref(:,1), ref(:,2), spline.tr, spline.control_points(1,:), spline.control_points(2,:));

step = 0.01;
spline.tr = spline.tr - step*dcost_dtr;
spline.control_points(1,:) = spline.control_points(1,:) - step*dcost_dPx;
spline.control_points(2,:) = spline.control_points(2,:) - step*dcost_dPy;

spline.tr = max(spline.tr, 1);
spline.tr = min(spline.tr, (spline.original_n-1+(spline.n-spline.original_n)/2));

% Evaluate Bspline with initial points
t = ((spline.n-spline.original_n)/2):0.02:(spline.original_n-1+(spline.n-spline.original_n)/2);
%t = 0:0.02:(spline.n-1);

pos = zeros(2, length(t));
for (i = 1:length(t))    
    pos(:,i) = spline.evaluate(t(i));
end

evaluation_points =  zeros(2, length(spline.tr));
for (i = 1:length(spline.tr))    
    evaluation_points(:,i) = spline.evaluate(spline.tr(i));
end

figure(1);
plot(spline.control_points(1,:), spline.control_points(2,:), 'o');
hold on;
plot(ref(:,1), ref(:,2), 'x');
plot(evaluation_points(1,:), evaluation_points(2,:), 'g^');
for (i = 1:length(spline.tr))   
    line([ref(i,1), evaluation_points(1,i)], [ref(i,2), evaluation_points(2,i)], 'Color', 'g', 'LineStyle', '--');
end
plot(pos(1,:), pos(2,:));
hold off;
axis equal;
xlim([0, 12]);
ylim([0, 10]);
if first
    drawnow;
    pause;
    first = false;
end
end


%% Try with fmincon
% Initial for fitting
spline = Bspline2.fit_reference(ref(:,1), ref(:,2), num_control_points, order);  
params0 = [spline.tr, spline.control_points(1,:), spline.control_points(2,:)]; % , 0:(size(spline.control_points, 2)-1)
tr1 = 1;
tr2 = tr1 + length(spline.tr) - 1;
Px1 = tr2 + 1;
Px2 = Px1 + size(spline.control_points, 2) - 1;
Py1 = Px2 + 1;
Py2 = Py1 + size(spline.control_points, 2) - 1;
tc1 = Py2 + 1;
tc2 = tc1 + size(spline.control_points, 2) - 1;

f_min = @(params) Bspline2.objective(order, ref(:,1), ref(:,2), params(tr1:tr2), params(Px1:Px2), params(Py1:Py2));  %params(tc1:tc2)
F_min = @(params) Bspline2.gradient2(order, ref(:,1), ref(:,2), params(tr1:tr2), params(Px1:Px2), params(Py1:Py2));
f_combined = @(params) deal(f_min(params), F_min(params));  % combine into function with multiple outputs

%opt.Jacobian='on';
%opt.Display='iter';
opt.Jacobian='romberg';
opt.Display='iter';
%[x1,~,fval_1,~,extra_arguments] = LevenbergMarquardt(f_min,params0,[],[],opt);


D1 = eye(length(spline.tr));
D2 = [zeros(size(D1,1),1), D1(:,1:end-1)];
D = D1 - D2;
A = [zeros(1,length(spline.tr)); D];
B = zeros(size(A,1), 1);
A(1,1) = -1;
B(1) = -1;
B(end) = (spline.original_n-1+(spline.n-spline.original_n)/2);

A = [A, zeros(size(A,1), length(params0)-size(A,2))];

Aeq = zeros(4, length(params0));
Beq = zeros(4, 1);
Aeq(1, Px1) = 1;
Aeq(1, Px1+1) = -1;
Aeq(2, Px2-1) = 1;
Aeq(2, Px2) = -1;
Aeq(3, Py1) = 1;
Aeq(3, Py1+1) = -1;
Aeq(4, Py2-1) = 1;
Aeq(4, Py2) = -1;

WithJacobian = false;
if (WithJacobian)
    options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'MaxFunctionEvaluations', 1e5, 'Display', 'iter', 'SpecifyObjectiveGradient', true); % run interior-point algorithm
    %x1 = fmincon(f_combined,params0,A,B,Aeq,Beq,[],[],[],options)
    x1 = fmincon(f_combined,params0,A,B,[],[],[],[],[],options)
else
    options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'MaxFunctionEvaluations', 1e5, 'Display', 'iter'); % run interior-point algorithm
    %x1 = fmincon(f_min,params0,A,B,Aeq,Beq,[],[],[],options)
    x1 = fmincon(f_min,params0,A,B,[],[],[],[],[],options)
end

%fmincon stopped because the size of the current step is less than
%the value of the step size tolerance and constraints are 
%satisfied to within the value of the constraint tolerance.
 


spline.tr = x1(tr1:tr2);
spline.control_points(1,:) = x1(Px1:Px2);
spline.control_points(2,:) = x1(Py1:Py2);

% Evaluate Bspline with initial points
t = ((spline.n-spline.original_n)/2):0.02:(spline.original_n-1+(spline.n-spline.original_n)/2);
%t = 0:0.02:(spline.n-1);

pos = zeros(2, length(t));
for (i = 1:length(t))    
    pos(:,i) = spline.evaluate(t(i));
end

evaluation_points =  zeros(2, length(spline.tr));
for (i = 1:length(spline.tr))    
    evaluation_points(:,i) = spline.evaluate(spline.tr(i));
end

figure(1);
plot(spline.control_points(1,:), spline.control_points(2,:), 'o');
hold on;
plot(ref(:,1), ref(:,2), 'x');
plot(evaluation_points(1,:), evaluation_points(2,:), 'g^');
for (i = 1:length(spline.tr))   
    line([ref(i,1), evaluation_points(1,i)], [ref(i,2), evaluation_points(2,i)], 'Color', 'g', 'LineStyle', '--');
end
plot(pos(1,:), pos(2,:));
hold off;
axis equal;
xlim([2, 8]);
ylim([0, 10]);


%%
figure(1);
clf;
xlim([0, 10]);
ylim([0, 10]);

ref = [];

while (true)
    current_ref = ginput(1);
    ref = [ref; current_ref];
    plot(ref(:,1), ref(:,2), 'x');
    xlim([0, 10]);
    ylim([0, 10]);

    if (size(ref,1) > 1)
        % Initial for fitting
        spline = Bspline2.fit_reference(ref(:,1), ref(:,2), length(ref), order);
        cost = spline.objective(order, ref(:,1), ref(:,2), spline.tr, spline.control_points(1,:), spline.control_points(2,:))
        
        % Evaluate Bspline with initial points
        t = ((spline.n-spline.original_n)/2):0.02:(spline.original_n-1+(spline.n-spline.original_n)/2);
        pos = zeros(2,length(t));        
        for (i = 1:length(t))    
            pos(:,i) = spline.evaluate(t(i));            
        end

        % figure(1);
        % plot(spline.control_points(1,:), spline.control_points(2,:), 'x');
        % xlim([0, 10]);
        % ylim([0, 10]);
        hold on;
        plot(pos(1,:), pos(2,:));
        hold off;
    end
end