function n = quat_rot(q,b)
%  n = quat_rot(q,b)
%  Rotate vector using quaternion
%
%   Input arguments:
%   q -  Attitude quaternion [1,4]
%   b -  Initial vector in body frame [3,1]
%
%   Output arguments:
%   n -  Rotated vector in navigation frame [3,1]

    Cbn = quat_dcm(q);
    n = Cbn*b;
end
