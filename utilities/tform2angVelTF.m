function B = tform2angVelTF(tform, sequence)
    if ~exist('sequence', 'var')
        % use the default ZYX axis sequence ...
        eul = WBM.utilities.tform2eul(tform);
    else
        eul = WBM.utilities.tform2eul(tform, sequence);
    end

    B = WBM.utilities.eul2angVelTF(eul);
end
