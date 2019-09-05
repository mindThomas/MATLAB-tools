function [pos, rotm] = frame2posRotm(vqT)
    if (size(vqT,1) ~= 7)
        error('frame2posRotm: %s', WBM.wbmErrorMsg.WRONG_VEC_DIM);
    end
    pos  = vqT(1:3,1);
    quat = vqT(4:7,1);

    % compute the orthonormal rotation matrix R ...
    rotm = WBM.utilities.quat2rotm(quat);
end
