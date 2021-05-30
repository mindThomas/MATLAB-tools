% 2D Rotation matrix
%
% R = ROT2(THETA) is a rotation matrix representing a rotation of THETA 
% radians about the z-axis.
%
% R = ROT2(THETA, 'deg') as above but THETA is in degrees.

function R = rot2(t, deg)
    if nargin > 1 && strcmp(deg, 'deg')
        t = t *pi/180;
    end
    
    ct = cos(t);
    st = sin(t);
    R = [
        ct  -st
        st   ct        
        ];
