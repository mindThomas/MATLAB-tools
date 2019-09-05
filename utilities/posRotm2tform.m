function tform = posRotm2tform(pos, rotm)
    if ( (size(pos,1) ~= 3) || (size(rotm,1) ~= 3) || (size(rotm,2) ~= 3) )
        error('posRotm2tform: %s', WBM.wbmErrorMsg.DIM_MISMATCH);
    end
    % create the homogeneous transformation matrix:
    tform = eye(4,4);
    tform(1:3,1:3) = rotm; % rotation
    tform(1:3,4)   = pos;  % translation
end
