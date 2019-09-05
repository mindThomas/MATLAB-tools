function axang = tform2axang(tform)
    if ( (size(tform,1) ~= 4) || (size(tform,2) ~= 4) )
        error('tform2axang: %s', WBM.wbmErrorMsg.WRONG_MAT_DIM);
    end
    % extract the rotation matrix and transform it ...
    R = tform(1:3,1:3);
    axang = WBM.utilities.rotm2axang(R);
end
