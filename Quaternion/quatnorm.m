function qout = quatnorm( q )
%  QUATNORM Calculate the norm of a quaternion.
%   N = QUATNORM( Q ) calculates the norm, N, for a given quaternion, Q.  Input
%   Q is an 4-by-M matrix containing M quaternions.  N returns a column vector
%   of M norms.  Each element of Q must be a real number.  Additionally, Q has
%   its scalar number as the first column.
%
%   Examples:
%
%   Determine the norm of q = [1 0 0 0]':
%      norm = quatnorm([1 0 0 0]')
%
%   See also QUATCONJ, QUATDIVIDE, QUATINV, QUATMOD, QUATMULTIPLY, 
%   QUATNORMALIZE, QUATROTATE.

if any(~isreal(q(:)))
    error('aero:quatnorm:isnotreal','Input elements are not real.');
end

if (size(q,1) ~= 4)
    error('aero:quatnorm:wrongdim','Input dimension is not 4-by-M.');
end

qout = sum(q.^2,1);
