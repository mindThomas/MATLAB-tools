function [R, t] = ICP2D(A, B, R0, t0)

    R2 = [R0, zeros(2,1); 0, 0, 1];
    t2 = [t0; 0];
    Btmp = R2(1:2,1:2)' * (B - t2(1:2));
    while (true)   
        [K, D] = dsearchn(Btmp', A');
        B_paired = Btmp(:,K);

        [Rtmp ttmp] = registerCensi(A, B_paired);
        %[Rtmp ttmp] = registerSVD([A; zeros(1,size(A,2))], [B_paired; zeros(1,size(B,2))]);
        R2 = R2 * Rtmp;
        t2 = t2 + ttmp;

        q = rotm2quat(Rtmp);
        if (2*acos(q(1)) < deg2rad(0.01) && norm(ttmp) < 0.00001) % stopping criteria
            break;
        end

        Btmp = R2(1:2,1:2)' * (B - t2(1:2));     
    end

    R = R2(1:2,1:2);
    t = t2(1:2);