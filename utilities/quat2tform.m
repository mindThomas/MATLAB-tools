function tform = quat2tform(quat)
    % build the homogeneous transformation matrix:
    tform = eye(4,4);
    tform(1:3,1:3) = WBM.utilities.quat2rotm(quat);
end
