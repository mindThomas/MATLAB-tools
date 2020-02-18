function yaw_error_corrected = HeadingErrorWrapping(yaw, yaw_ref)
    err = yaw_ref - yaw;
    if (err >= 0 && err < pi)
        yaw_error_corrected = mod(err + pi, pi);
    else
        yaw_error_corrected = mod(err + pi, pi) - pi;
    end
end
