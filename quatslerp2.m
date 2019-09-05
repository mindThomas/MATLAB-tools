% Quaternion slerp for rotation quaternions
% t is the progress parameter from 0 to 1
% https://ww2.mathworks.cn/help/fusion/examples/lowpass-filter-orientation-using-quaternion-slerp.html
function q_out = quatslerp2(q0, q1, t)
    % Compute angle between the two rotations
    %q_diff = quatmultiply(quatconj(q0), q1);
    %angle_diff = 2 * acos(q_diff(1));
    % This is however also the same as
    angle_diff = 2 * acos( dot(q0, q1) );

    % Compute the slerped quaternion
    q_out = q0 * sin( 1/2 * angle_diff * (1-t) ) / sin(1/2 * angle_diff) ...
          + q1 * sin ( 1/2 * angle_diff * t ) / sin(1/2 * angle_diff);
end
