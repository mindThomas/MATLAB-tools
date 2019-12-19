function [Euler_, bw_, Angles_] = test_ahrs_dcm_data(dWb, Fb, Mb, fn, mn, dt)
%  [Euler_, bw_, Angles_] = test_ahrs_dcm_data(dWb, Fb, Mb, q_ref, q_init, dt)
%  Runs the sumulation for the DCM-based AHRS algorithm with sensors data
%
%   Input arguments:
%   dWb -  Gyroscopes measurements (integral of the angular rate over 
%   the computer cycle)
%   Fb  - Accelerometers measurements
%   Mb  - Magnetometer measurements
%   fn - Gravity vector in navigation frame
%   mn - Magnetic field vector in navigation frame 
%   dt - time interval between measurements (computer cycle)
%
%   Output arguments:
%   Euler_  - Estimated attitude angles (AHRS)
%   bw_     - Estimated gyroscopes biases
%   Euler_  - Estimated attitude angles (Accelerometers only)

%%
close all; clc;

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

%% INS Axes
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
Cbn = eye(3);

%Gyroscope bias estimate
bw = [0; 0; 0];

%Initial Covariance matrix
P=zeros(6,6);
P(1:3,1:3)=diag([1e-2, 1e-2, 1e-2]); 
P(4:6,4:6)=diag([1e-6, 1e-6, 1e-6]); 

%% Logs
Nsim = size(dWb,1);  %number of simulation steps
Euler_ = zeros(Nsim,3);
bw_ = zeros(Nsim,3);
Angles_ = zeros(Nsim,2);

%% Main loop
for i=1:Nsim
    
    
    %% Inertial Sensors readings with noises and biases
    dwb = dWb(i,:)';
    fb  = Fb(i,:)'./norm(Fb(i,:));
    mb  = Mb(i,:)'./norm(Mb(i,:));
    
    %% Call AHRS Algorithm
    [Cbn, P, bw] = ahrs_dcm(Cbn, P, bw, dwb, fb, mb, fn, mn, dt);
    
    %% collect data logs
    bw_(i,:) = bw*dt;
    [Euler_(i,1), Euler_(i,2), Euler_(i,3)] = dcm_angle(Cbn');
    
    %% Accelerometer-based angles
    Angles_(i,1) =  atan2(fb(1),sqrt(fb(2)^2+fb(3)^2));
    Angles_(i,2) = -atan2(fb(2),sqrt(fb(1)^2+fb(3)^2));
    
    %% Animations
    if(mod(i,animation_rate) == 0)
        
        %% INS
        %axes
        %body and body axes
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
plot(t,Angles_*180/pi);
plot(t,Euler_*180/pi,'Linewidth',1);
legend('\theta_{accl}','\gamma_{accl}','\psi_{ahrs}','\theta_{ahrs}','\gamma_{ahrs}');
ylabel('Angles, deg');
xlabel('Time,sec');

figure('Name','Gyroscopes biases estimates');
hold on;
grid on;
plot(t,bw_*dt);
legend('\delta\omega_x','\delta\omega_y','\delta\omega_z');
ylabel('Gyro Bias, rad/sec');
xlabel('Time,sec');

figure('Name','Attitude Comparison');
subplot(2,1,1);
hold on;
grid on;
plot(t,Angles_(:,1)*180/pi,'r--');
plot(t,Euler_(:,2)*180/pi,'r-','Linewidth',1);
plot(t,Angles_(:,2)*180/pi,'b--');
plot(t,Euler_(:,3)*180/pi,'b-','Linewidth',1);
legend('r^{accl}_y','r^{ahrs}_y','r^{accl}_x','r^{ahrs}_x');
ylabel('Angles, deg');
xlabel('Time,sec');
axis tight;
subplot(2,1,2);
hold on;
grid on;
plot(t,sqrt(Fb(:,1).^2+Fb(:,2).^2+Fb(:,3).^2));
ylabel('Acceleration norm, m/sec^2');
xlabel('Time,sec');
axis tight;

end

