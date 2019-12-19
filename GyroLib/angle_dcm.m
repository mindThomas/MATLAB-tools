function Cbn = angle_dcm(rz, ry, rx)
%  Cbn = angle_dcm(rz, ry, rx)
%  Transforms Euler angles to  direction cosine matrix
%
%   Input arguments:
%   rz, ry, rx - Euler angles around Z, Y and X axes correspondingly
%
%   Output arguments:
%   Cbn -  Direction cosine matrix [3,3]

sz = sin(rz);
cz = cos(rz);
sy = sin(ry);
cy = cos(ry);
sx = sin(rx);
cx = cos(rx);

Cx = [1 0 0; 0 cx sx; 0 -sx cx];
Cy = [cy 0 -sy; 0 1 0; sy 0 cy];
Cz = [cz sz 0; -sz cz 0; 0 0 1];

Cbn = Cx*Cy*Cz;
