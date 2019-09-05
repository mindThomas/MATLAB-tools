% Take natural logarithm of a quaternion
% https://math.stackexchange.com/a/939288
function q_out = quatlogn(q)
    q0 = q(1);
    q_vec = q(2:4);
    norm_q = norm(q_vec);

    if (norm_q > eps)
        q_out = [0; acos(q0) * q_vec / norm_q];
    else
        q_out = [0; 0; 0; 0]; % return error/undefined quaternion
    end
end