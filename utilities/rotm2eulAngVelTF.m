function [eul, B] = rotm2eulAngVelTF(rotm, sequence)
    if ~exist('sequence', 'var')
        % use the default sequence ...
        sequence = 'ZYX';
    end
    eul = WBM.utilities.rotm2eul(rotm, sequence);
    B   = WBM.utilities.eul2angVelTF(eul, sequence);
end
