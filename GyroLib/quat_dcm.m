function Cbn = quat_dcm(q)
%  Cbn = quat_dcm(q)
%  Transform quaternion to direction cosine matrix
%
%   Input arguments:
%   q -  Attitude quaternion [1,4]
%
%   Output arguments:
%   Cbn - direction cosine matrix

q0 = q(1);
q1 = q(2);
q2 = q(3);
q3 = q(4);

Cbn_1_1 = q0^2+q1^2-q2^2-q3^2;
Cbn_1_2 = 2*(q1*q2+q0*q3);
Cbn_1_3 = 2*(q1*q3-q0*q2);
Cbn_2_1 = 2*(q1*q2-q0*q3);
Cbn_2_2 = q0^2-q1^2+q2^2-q3^2;
Cbn_2_3 = 2*(q2*q3+q0*q1);
Cbn_3_1 = 2*(q1*q3+q0*q2);
Cbn_3_2 = 2*(q2*q3-q0*q1);
Cbn_3_3 = q0^2-q1^2-q2^2+q3^2;

Cbn = [
    Cbn_1_1 Cbn_1_2 Cbn_1_3
    Cbn_2_1 Cbn_2_2 Cbn_2_3
    Cbn_3_1 Cbn_3_2 Cbn_3_3
    ];

end
