function eul = tform2eul(tform, sequence)
    if ( (size(tform,1) ~= 4) || (size(tform,2) ~= 4) )
        error('tform2eul: %s', WBM.wbmErrorMsg.WRONG_MAT_DIM);
    end
    % extract the rotation matrix and transform it ...
    R = tform(1:3,1:3);
    if ~exist('sequence', 'var')
        % use the default axis sequence ZYX ...
        eul = WBM.utilities.rotm2eul(R);
        return
    end
    % else ...
    eul = WBM.utilities.rotm2eul(R, sequence);
end
