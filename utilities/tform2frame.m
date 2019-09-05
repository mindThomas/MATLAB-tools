function vqT = tform2frame(tform)
    if ( (size(tform,1) ~= 4) || (size(tform,2) ~= 4) )
        error('tform2frame: %s', WBM.wbmErrorMsg.WRONG_MAT_DIM);
    end
    vqT = zeros(7,1);
    R   = tform(1:3,1:3); % extract the rotation matrix ...

    vqT(1:3,1) = tform(1:3,4);               % translation
    vqT(4:7,1) = WBM.utilities.rotm2quat(R); % orientation
end
