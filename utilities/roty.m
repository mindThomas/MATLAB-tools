function R_y = roty(ang)
    R_y = eye(3,3);

    sa = sin(ang);
    ca = cos(ang);

    R_y(1,1) =  ca;
    R_y(1,3) =  sa;
    R_y(3,1) = -sa;
    R_y(3,3) =  ca;
end
