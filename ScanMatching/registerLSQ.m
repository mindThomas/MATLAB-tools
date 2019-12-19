function [R,t] = registerLSQ(P1, P2)   
    % Find rotation and translation to 
    % match points in P2 ~ [2xm]
    % to the points in P1 ~ [2xn]

    % Make sure that P2 is the smallest   
    if (size(P2, 2) > size(P1, 2))
        P2tmp = P2;
        P2 = P1;
        P1 = P2tmp;
    end    
    
    function [v] = fmin(X)
        rot = [cos(X(3)) -sin(X(3)); sin(X(3)) cos(X(3))];
        
        % Apply transform
        P2_prime = rot * P2 + [X(1); X(2)];

        % Compute associations/correspondence pairs
        % Compute a closest index in P1 for each point in P2
        [K, D] = dsearchn(P1', P2_prime');
        
        % Construct output vector of vectorial difference
        v = zeros(2*length(K),1);        
        for k = 0:(length(K)-1)
            v(2*k+1:2*k+2) = P1(:, K(k+1)) - P2(:, k+1);
        end

        v
        
    end

    x0 = 0;
    y0 = 0;
    yaw0 = 0;

    [T res] = lsqnonlin(@fmin, [x0 y0 yaw0]);
    T = T - [x0 y0 yaw0];
    R = rotz(T(3));
    t = [T(1:2) 0];
end