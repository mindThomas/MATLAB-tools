function [rz, ry, rx] = quat_angle(q)
%  [rz, ry, rx] = quat_angle(q)
%  Transforms quaternion to Euler angles
%
%   Input arguments:
%   q -  Attitude quaternion [1,4]
%
%   Output arguments:
%   rz, ry, rx - Euler angles around Z, Y and X axes correspondingly

Cbn = quat_dcm(q);
[rz, ry, rx] = dcm_angle(Cbn);

end