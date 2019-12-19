function q = dcm_quat(Cbn)
%  Cbn = quat_dcm(q)
%  Transform  direction cosine matrix to quaternion
%
%   Input arguments:
%   Cbn - direction cosine matrix
%
%   Output arguments:
%   q -  Attitude quaternion [1,4]


q0 = 1/2*(1+Cbn(1,1)+Cbn(2,2)+Cbn(3,3))^(1/2);
q1 = (1/(4*q0))*(Cbn(2,3)-Cbn(3,2));
q2 = (1/(4*q0))*(Cbn(3,1)-Cbn(1,3));
q3 = (1/(4*q0))*(Cbn(1,2)-Cbn(2,1));

q = [q0, q1, q2, q3];

end
