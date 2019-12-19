function [Euler_, bw_] = test_ahrs_quat(dWb, Fb, Mb, q_ref, q_init, fn, mn, dt)
%  [Euler_, bw_] = test_ahrs_quat(dWb, Fb, Mb, q_ref, q_init, dt)
%  Runs the sumulation for the quaternion-based AHRS algorithm
%
%   Input arguments:
%   dWb -  Gyroscopes measurements (integral of the angular rate over 
%   the computer cycle)
%   Fb  - Accelerometers measurements
%   Mb  - Magnetometer measurements
%   q_ref  - Reference attitude quaternion readings
%   q_init - Initial reference quaternion
%   fn - Gravity vector in navigation frame
%   mn - Magnetic field vector in navigation frame 
%   dt - time interval between measurements (computer cycle)
%
%   Output arguments:
%   Euler_  - Estimated attitude angles
%   bw_     - Estimated gyroscopes biases

%%
close all; clc;

%% randomize
rng('shuffle');

%% figure
figure;
view(3);
axis equal;
hold on; grid on;
set(gca,'Xlim',[-2 2]);
set(gca,'Ylim',[-2 2]);
set(gca,'Zlim',[-2 2]);
set(gca,'YDir','reverse');
set(gca,'ZDir','reverse');
set(gcf,'renderer','opengl');
animation_rate = 30;

%% Target axes
at1 = [0; 0; 2.5]*4;
at2 = [0; 2.5; 0]*4;
at3 = [2.5; 0; 0]*4;
st =  [0; 0; 0]*3;
pat1 = plot3([st(1) at1(1)],[st(2) at1(2)], [st(3) at1(3)],'b--','linewidth',0.1);
pat2 = plot3([st(1) at2(1)],[st(2) at2(2)], [st(3) at2(3)],'g--','linewidth',0.1);
pat3 = plot3([st(1) at3(1)],[st(2) at3(2)], [st(3) at3(3)],'r--','linewidth',0.1);

%% AHRS Axes
ai1 = [0; 0; 1.5;];
ai2 = [0; 1.5; 0;];
ai3 = [1.5; 0; 0;];
sz =  [0; 0; 0;];
pai1 = ...
    plot3([sz(1) ai1(1)],[sz(2) ai1(2)], [sz(3) ai1(3)],'b-','linewidth',3);
pai2 = ...
    plot3([sz(1) ai2(1)],[sz(2) ai2(2)], [sz(3) ai2(3)],'g-','linewidth',3);
pai3 = ...
    plot3([sz(1) ai3(1)],[sz(2) ai3(2)], [sz(3) ai3(3)],'r-','linewidth',3);

%% Navigation frame Axes
an1 = [0; 0; 10;];
an2 = [0; 10; 0;];
an3 = [10; 0; 0;];
sz =  [0; 0; 0;];
plot3([sz(1) an1(1)],[sz(2) an1(2)], [sz(3) an1(3)],'b-','linewidth',0.1);
plot3([sz(1) an2(1)],[sz(2) an2(2)], [sz(3) an2(3)],'g-','linewidth',0.1);
plot3([sz(1) an3(1)],[sz(2) an3(2)], [sz(3) an3(3)],'r-','linewidth',0.1);
an1 = [0; 0; -10;];
an2 = [0; -10; 0;];
an3 = [-10; 0; 0;];
sz =  [0; 0; 0;];
plot3([sz(1) an1(1)],[sz(2) an1(2)], [sz(3) an1(3)],'b-','linewidth',0.1);
plot3([sz(1) an2(1)],[sz(2) an2(2)], [sz(3) an2(3)],'g-','linewidth',0.1);
plot3([sz(1) an3(1)],[sz(2) an3(2)], [sz(3) an3(3)],'r-','linewidth',0.1);

%% Cube patch for AHRS attitude animation
Vert = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1];
Vert(:,1) = (Vert(:,1)-0.5);
Vert(:,2) = (Vert(:,2)-0.5);
Vert(:,3) = (Vert(:,3)-0.5);
Faces = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
ptch.Vertices = Vert;
ptch.Faces = Faces;
ptch.FaceVertexCData = copper(6);
ptch.FaceColor = 'flat';
patch_handle = patch(ptch);

%% Text
th = text(-7,1,'Time ');
set(th,'FontSize',9,'Color','b','FontWeight','bold');

%% Initial values
%Attitude
q = q_init;
%Gyroscope bias estimate
bw = [0; 0; 0];

%Introduce initial attitude errors
att_err = 1e-1;
err_x = randn*sqrt(att_err);
err_y = randn*sqrt(att_err);
err_z = randn*sqrt(att_err);
qe = [1, err_x/2, err_y/2, err_z/2];
qe = qe/sqrt(qe*qe');
q = quat_mult( q, qe );

%Initial Covariance matrix
P=zeros(6,6);
P(1:3,1:3)=diag([1e-2, 1e-2, 1e-2]); %attitude errors in radian
P(4:6,4:6)=diag([1e-6, 1e-6, 1e-6]); %gyro bias errors in rad/s


%% Logs
Nsim = size(dWb,1);  %number of simulation steps
Euler_ = zeros(Nsim,3);
bw_ = zeros(Nsim,3);
err_ = zeros(Nsim,3);

%% Main loop
for i=1:Nsim
    
    
    %% Inertial Sensors readings with noises and biases
    dwb = dWb(i,:)';
    fb = Fb(i,:)';
    mb = Mb(i,:)';
    
    %% Call AHRS Algoryithm
    [q, P, bw] = ahrs_quat(q, P, bw, dwb, fb, mb, fn, mn, dt);
    
    %% Collect data logs
    bw_(i,:) = bw*dt;
    [Euler_(i,1), Euler_(i,2), Euler_(i,3)] = quat_angle(quat_conj(q));
    
    %Angles errors
    qerr = quat_mult(quat_conj(q),q_ref(i,:));
    [err_(i,1), err_(i,2), err_(i,3)] = quat_angle(qerr); 
    
    %% Animation
    if(mod(i,animation_rate) == 0)
        
        %Reference axes
        Cbnref = quat_dcm(q_ref(i,:))';
        a1_ = Cbnref'*at1;
        a2_ = Cbnref'*at2;
        a3_ = Cbnref'*at3;
        s_ = [0,0,0];
        set(pat1,'Xdata',[s_(1) a1_(1)]);
        set(pat1,'Ydata',[s_(2) a1_(2)]);
        set(pat1,'Zdata',[s_(3) a1_(3)]);
        set(pat2,'Xdata',[s_(1) a2_(1)]);
        set(pat2,'Ydata',[s_(2) a2_(2)]);
        set(pat2,'Zdata',[s_(3) a2_(3)]);
        set(pat3,'Xdata',[s_(1) a3_(1)]);
        set(pat3,'Ydata',[s_(2) a3_(2)]);
        set(pat3,'Zdata',[s_(3) a3_(3)]);
        
        %AHRS axes
        Cbn = quat_dcm(q);
        a1_ = Cbn*ai1;
        a2_ = Cbn*ai2;
        a3_ = Cbn*ai3;
        s_ = [0, 0, 0];
        set(pai1,'Xdata',[s_(1) a1_(1)]);
        set(pai1,'Ydata',[s_(2) a1_(2)]);
        set(pai1,'Zdata',[s_(3) a1_(3)]);
        set(pai2,'Xdata',[s_(1) a2_(1)]);
        set(pai2,'Ydata',[s_(2) a2_(2)]);
        set(pai2,'Zdata',[s_(3) a2_(3)]);
        set(pai3,'Xdata',[s_(1) a3_(1)]);
        set(pai3,'Ydata',[s_(2) a3_(2)]);
        set(pai3,'Zdata',[s_(3) a3_(3)]);
        
        %Cube
        Vert_ = Vert;
        for j=1:size(Vert,1)
            Vert_(j,:) = (Vert(j,:)*Cbn');
        end
        set(patch_handle,'Vertices',Vert_);
        
        set(th,'String',...
            sprintf('tm %2.1f\nps   % 7.5f\ntt     % 7.5f\ngm  % 7.5f',...
            i*dt,Euler_(i,1)*180/pi,Euler_(i,2)*180/pi,Euler_(i,3)*180/pi));
        
        drawnow;
    end
end

%% Plot Results
t = (1:Nsim)*dt;
figure('Name','Attitude');
hold on;
grid on;
plot(t,Euler_*180/pi);
legend('\psi','\theta','\gamma');
ylabel('Aangles, deg');
xlabel('Time,sec');
figure('Name','Inertial Sensors Errors estimates');
hold on;
grid on;
plot(t,bw_*dt);
legend('\delta\Theta_x','\delta\Theta_y','\delta\Theta_z');
ylabel('Gyro Bias, rad/sec');
xlabel('Time,sec');
figure('Name','Angles Errors');
hold on;
grid on;
plot(t,err_*180/pi);
legend('\delta\psi','\delta\theta','\delta\gamma');
ylabel('Angles errors, deg');
xlabel('Time,sec');

end


