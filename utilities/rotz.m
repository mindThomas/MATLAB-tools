function R_z = rotz(ang)
    R_z = eye(3,3);

    sa = sin(ang);
    ca = cos(ang);

    R_z(1,1) =  ca;
    R_z(1,2) = -sa;
    R_z(2,1) =  sa;
    R_z(2,2) =  ca;
end
