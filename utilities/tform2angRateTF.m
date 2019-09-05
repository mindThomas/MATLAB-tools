function B_inv = tform2angRateTF(tform, sequence)
    if ~exist('sequence', 'var')
        % use the default ZYX axis sequence ...
        eul = WBM.utilities.tform2eul(tform);
    else
        eul = WBM.utilities.tform2eul(tform, sequence);
    end

    B_inv = WBM.utilities.eul2angRateTF(eul);
end
