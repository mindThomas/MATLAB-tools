% Quaternion exponential
% See appendix H.5 in Kugle - Modelling and Control of a Ball-balancing Robot
function q_out = quatexp(q)
    q0 = q(1);
    q_vec = q(2:4);
    norm_q = norm(q_vec);

    if (norm_q > eps)
        q_out = exp(q0) * [cos(norm_q); sin(norm_q) * q_vec / norm_q];
    else
        q_out = [1; 0; 0; 0]; % return unit quaternion
    end
end