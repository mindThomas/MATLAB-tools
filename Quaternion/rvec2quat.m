function q = rvec2quat(angle, rvec)
        rvec_ = [rvec(1); rvec(2); rvec(3)];
        q = [cos(angle/2); rvec_ * sin(angle/2)];
end
