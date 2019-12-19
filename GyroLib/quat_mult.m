function q = quat_mult(lam, mu)
%  q = quat_mult(lam, mu)
%  Multiplies quaternions  q = lam*mu
%
%   Input arguments:
%   lam -  Attitude quaternion [1,4]
%   mu -  Attitude quaternion [1,4]
%
%   Output arguments:
%   q -   Attitude quaternion [1,4]

l0 = lam(1);
l1 = lam(2);
l2 = lam(3);
l3 = lam(4);

m = [
    l0 -l1 -l2 -l3
    l1  l0 -l3  l2
    l2  l3  l0 -l1
    l3 -l2  l1  l0
    ];
q = (m*mu')';
    
end
