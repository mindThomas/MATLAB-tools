% Quaternion slerp for rotation quaternions
% t is the progress parameter from 0 to 1
% https://ww2.mathworks.cn/help/fusion/examples/lowpass-filter-orientation-using-quaternion-slerp.html
function q_out = quatslerp(q0, q1, t)
    % Compute difference quaternion
    q_diff = quatmultiply(quatconj(q0), q1);

    % Compute delta quaternion based on progress parameter and quaternion power
    q_delta = quatpow(q_diff, t);

    % Add delta quaternion to start quaternion
    q_out = quatmultiply(q0, q_delta);
end
