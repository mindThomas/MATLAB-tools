function [R, t] = ScanICP(A, B, R, t)
    R2 = [R, zeros(2,1); 0,0,1];
    t2 = [t; 0];

    Btmp = R2(1:2,1:2)' * (B - t2(1:2));
    while (true)       
        [K, D] = dsearchn(Btmp', A');
        B_paired = Btmp(:,K);

        excludedIdx = find(D > 1); % only include points which are no further apart than 1 meter

        Aprime = A;
        Bprime = B_paired;
        Aprime(:,excludedIdx) = [];
        Bprime(:,excludedIdx) = [];

        [Rtmp ttmp] = registerCensi(Aprime, Bprime);
        %[Rtmp ttmp] = registerSVD([A; zeros(1,size(A,2))], [B_paired; zeros(1,size(B,2))])
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
        for (j = 1:size(Aprime,2))
            line([Aprime(1,j), Bprime(1,j)], [Aprime(2,j), Bprime(2,j)]);
        end
        hold off;
        axis equal;    
        pause(0.1);
    end
    
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
    
    R = R2(1:2,1:2);
    t = t2(1:2);
end