function test_spin_cone( dWb, Fb, Mb, fn, mn, q_ref, q_init, dt )

close all;

%% figure
figure;
view(3);
axis equal;
hold on; grid on;
set(gca,'Xlim',[-10 10]);
set(gca,'Ylim',[-10 10]);
set(gca,'Zlim',[-10 10]);
set(gca,'YDir','reverse');
set(gca,'ZDir','reverse');
set(gcf,'renderer','opengl');
animation_rate = 1/dt/5;

%% Target axes
at1 = [0; 0; 2.5]*4;
at2 = [0; 2.5; 0]*4;
at3 = [2.5; 0; 0]*4;
st =  [0; 0; 0]*3;
pat1 = plot3([st(1) at1(1)],[st(2) at1(2)], [st(3) at1(3)],'b--','linewidth',0.1);
pat2 = plot3([st(1) at2(1)],[st(2) at2(2)], [st(3) at2(3)],'g--','linewidth',0.1);
pat3 = plot3([st(1) at3(1)],[st(2) at3(2)], [st(3) at3(3)],'r--','linewidth',0.1);

%% INS Axes
ai1 = [0; 0; 1.5;]*3;
ai2 = [0; 1.5; 0;]*3;
ai3 = [1.5; 0; 0;]*3;
sz =  [0; 0; 0;]*3;
pai1 = plot3([sz(1) ai1(1)],[sz(2) ai1(2)], [sz(3) ai1(3)],'b-','linewidth',3);
pai2 = plot3([sz(1) ai2(1)],[sz(2) ai2(2)], [sz(3) ai2(3)],'g-','linewidth',3);
pai3 = plot3([sz(1) ai3(1)],[sz(2) ai3(2)], [sz(3) ai3(3)],'r-','linewidth',3);

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

%% Cube patch for INS attitude animation
Vert = [ ...
    0 0 0 ;
    1 0 0 ;
    1 1 0 ;
    0 1 0 ;
    0 0 1 ;
    1 0 1 ;
    1 1 1 ;
    0 1 1];
Vert(:,1) = (Vert(:,1)-0.5)*3;
Vert(:,2) = (Vert(:,2)-0.5)*3;
Vert(:,3) = (Vert(:,3)-0.5)*3;
Faces = [ ...
    1 2 6 5 ;
    2 3 7 6 ;
    3 4 8 7 ;
    4 1 5 8 ;
    1 2 3 4 ;
    5 6 7 8 ];
ptch.Vertices = Vert;
ptch.Faces = Faces;
ptch.FaceVertexCData = ones(6,3)*0.7;
ptch.FaceColor = 'flat';
patch_handle = patch(ptch);

th = text(-2,20,20,'Time ');
set(th,'FontSize',12,'Color','b');

%% Initial values
%Initial alignment
q  = q_init;

%% Logs
Nsim = size(q_ref,1); %number of simulation steps
Angles  = zeros(Nsim,3);
Euler = zeros(Nsim,3);
dEuler = zeros(Nsim,3);
Quest = zeros(Nsim,3);
dQuest = zeros(Nsim,3);
Triad = zeros(Nsim,3);
dTriad = zeros(Nsim,3);

%% Main loop
for i=1:Nsim
    
    %% Gyroscope
    gamma1 = dWb(i,1);
    gamma2 = dWb(i,2);
    gamma3 = dWb(i,3);
    gamma = norm(dWb(i,:));
    lambda0 =cos(gamma/2);
    lambda1 = gamma1*sin(gamma/2)/gamma;
    lambda2 = gamma2*sin(gamma/2)/gamma;
    lambda3 = gamma3*sin(gamma/2)/gamma;
    if (gamma~=0)
        lambda = [lambda0, lambda1, lambda2, lambda3];
    else
        lambda = [1, 0, 0, 0];
    end
    q = quat_mult(quat_conj(lambda),q);
    q = q/sqrt(q*q');
    Cbn_gyro = quat_dcm(q);
    
    
    %% TRIAD
    Cbn_triad = TRIAD(Fb(i,:)', Mb(i,:)', fn, mn);
    [Triad(i,1), Triad(i,2), Triad(i,3)] = ...
        dcm_angle(Cbn_triad');
    
    %% QUEST
    wf = 1.0;
    wm = 1.0;
    quest = QUEST(Fb(i,:)', Mb(i,:)', fn, mn, wf, wm);
    [Quest(i,1), Quest(i,2), Quest(i,3)] = ...
        quat_angle(quat_conj(quest));
    
    %% Collect data logs
    %Attitude angles
    [Euler(i,1), Euler(i,2), Euler(i,3)] = quat_angle(quat_conj(q));
    [Angles(i,1), Angles(i,2), Angles(i,3)] = ...
        quat_angle(quatconj(q_ref(i,:)));
    
    %Gyroscopes angles errors
    Cref = quat_dcm(q_ref(i,:));
    Cimu = quat_dcm(q);
    Cerr = Cref*Cimu';
    [dEuler(i,1), dEuler(i,2), dEuler(i,3)] = dcm_angle(Cerr);
    
    %TRIAD angles errors
    Cref = quat_dcm(q_ref(i,:));
    Cerr = Cref*Cbn_triad';
    [dTriad(i,1), dTriad(i,2), dTriad(i,3)] = dcm_angle(Cerr);
    
    %QUEST angles errors
    Cref = quat_dcm(q_ref(i,:));
    Cimu = quat_dcm(quest);
    Cerr = Cref*Cimu';
    [dQuest(i,1), dQuest(i,2), dQuest(i,3)] = dcm_angle(Cerr);
    
    %% Animations
    if(mod(i,animation_rate) == 0)
        
        %% Reference
        %trajectory
        %body and body axes
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
        
        %% Estimated
        %axes
        %body and body axes
        a1_ = Cbn_gyro*ai1;
        a2_ = Cbn_gyro*ai2;
        a3_ = Cbn_gyro*ai3;
        s_ = [0,0,0];
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
            Vert_(j,:) = (Vert(j,:)*Cbn_gyro');
        end
        set(patch_handle,'Vertices',Vert_);
        
        drawnow;
    end
end

%% Plot Results
t = (1:Nsim)*dt;
figure('Name','Attitude Gyroscope');
subplot(2,1,1);hold on; grid on;
plot(t,Euler);
plot(t,Angles,'--');
legend('Z_{est}','Y_{est}','X_{est}','Z_{ref}','Y_{ref}','X_{ref}');
ylabel('Angles, rad');
xlabel('Time,sec');
subplot(2,1,2);hold on; grid on;
plot(t,dEuler*180/pi);
ylabel('Angles error, deg');
xlabel('Time,sec');

figure('Name','Attitude QUEST');
subplot(2,1,1);hold on; grid on;
plot(t,Quest);
plot(t,Angles,'--');
legend('Z_{est}','Y_{est}','X_{est}','Z_{ref}','Y_{ref}','X_{ref}');
ylabel('Angles, rad');
xlabel('Time,sec');
subplot(2,1,2);hold on; grid on;
plot(t,dQuest*180/pi);
ylabel('Angles error, deg');
xlabel('Time,sec');

figure('Name','Attitude TRIAD');
subplot(2,1,1);hold on; grid on;
plot(t,Triad);
plot(t,Angles,'--');
legend('Z_{est}','Y_{est}','X_{est}','Z_{ref}','Y_{ref}','X_{ref}');
ylabel('Angles, rad');
xlabel('Time,sec');
subplot(2,1,2);hold on; grid on;
plot(t,dTriad*180/pi);
ylabel('Angles error, deg');
xlabel('Time,sec');
