A = 20*rand(2, 100) - 10;

R = rot2(deg2rad(10));
t = [0.5; 0.7];

B = R*A + t + 0.1*randn(2, 100);

B = B(:, randperm(size(B, 2)));

R2 = eye(3);
t2 = [0;0;0];
Btmp = R2(1:2,1:2)' * (B - t2(1:2));
while (true)   
    [K, D] = dsearchn(Btmp', A');
    B_paired = Btmp(:,K);
    
    %[Rtmp ttmp] = registerCensi(A, B_paired);
    [Rtmp ttmp] = registerSVD([A; zeros(1,size(A,2))], [B_paired; zeros(1,size(B,2))])
    R2 = R2 * Rtmp;
    t2 = t2 + ttmp;
    
    q = rotm2quat(Rtmp);
    if (2*acos(q(1)) < deg2rad(0.01) && norm(ttmp) < 0.00001) % stopping criteria
        break;
    end
    
    Btmp = R2(1:2,1:2)' * (B - t2(1:2));
    
    figure(1);    
    plot(A(1,:), A(2,:), 'r.');
    hold on;
    plot(B_paired(1,:), B_paired(2,:), 'b.');
    for (j = 1:size(A,2))
        line([A(1,j), B_paired(1,j)], [A(2,j), B_paired(2,j)]);
    end
    hold off;
    axis equal;        
end

rad2deg(rotm2eul(R2))*[1,0,0]'
C = R2(1:2,1:2)' * (B - t2(1:2));
t2(1:2)

figure(2);
subplot(1,2,1);
plot(A(1,:), A(2,:), 'r.');
hold on;
plot(B(1,:), B(2,:), 'b.');
hold off;
axis equal;

subplot(1,2,2);
plot(A(1,:), A(2,:), 'r.');
hold on;
plot(C(1,:), C(2,:), 'b.');
hold off;
axis equal;