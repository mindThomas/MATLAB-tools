function qconj = quat_conj(q)
%  q = quat_conj(q)
%  Conjugate quaternion
%
%   Input arguments:
%   q -  Attitude quaternion [1,4]
%
%   Output arguments:
%   qconj -  Conjugate attitude quaternion [1,4]

%% Compute conjugate quaternion
q0 = q(1);
q1 = q(2);
q2 = q(3);
q3 = q(4);

qconj = [q0, -q1, -q2, -q3];
end
