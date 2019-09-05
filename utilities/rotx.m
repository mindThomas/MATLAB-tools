function R_x = rotx(ang)
    R_x = eye(3,3);

    sa = sin(ang);
    ca = cos(ang);

    R_x(2,2) =  ca;
    R_x(2,3) = -sa;
    R_x(3,2) =  sa;
    R_x(3,3) =  ca;
end
