function out = quatmul(q, p)
    % for q o p = Phi(q) * p
    Phi_q = [q(1) -q(2) -q(3) -q(4);     % for q o p = Phi(q) * p
             q(2) q(1)  -q(4) q(3);
             q(3) q(4)  q(1)  -q(2);
             q(4) -q(3) q(2)  q(1)]
    
    out = Phi_q * p;