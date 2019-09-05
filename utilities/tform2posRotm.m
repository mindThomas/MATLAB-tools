function [pos, rotm] = tform2posRotm(tform)
    if ( (size(tform,1) ~= 4) || (size(tform,2) ~= 4) )
        error('tform2posRotm: %s', WBM.wbmErrorMsg.WRONG_MAT_DIM);
    end
    % extract the translation vector and the rotation matrix ...
    pos  = tform(1:3,4);
    rotm = tform(1:3,1:3);
end
