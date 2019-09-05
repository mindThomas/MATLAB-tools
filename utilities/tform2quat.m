function quat = tform2quat(tform)
    if ( (size(tform,1) ~= 4) || (size(tform,2) ~= 4) )
        error('tform2quat: %s', WBM.wbmErrorMsg.WRONG_MAT_DIM);
    end
    % extract the rotation matrix and transform it ...
    R = tform(1:3,1:3);
    quat = WBM.utilities.rotm2quat(R);
end
