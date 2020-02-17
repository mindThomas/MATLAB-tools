A = 20*rand(2, 100) - 10;

R = rot2(deg2rad(35));
t = [1; 2];

B = R*A + t + 0.5*randn(2, 100);

[R2 t2] = registerCensi(A, B)
[R3 t3] = registerSVD([A; zeros(1,size(A,2))], [B; zeros(1,size(B,2))])
%[R3 t3] = registerLSQ(A, B)
rad2deg(rotm2eul(R2))*[1,0,0]'
C = R2(1:2,1:2)' * (B - t2(1:2));


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