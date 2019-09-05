function [pos, eul] = frame2posEul(vqT, sequence)
     if ~exist('sequence', 'var')
        % use the default sequence ...
        sequence = 'ZYX';
    end
    [pos, R] = WBM.utilities.frame2posRotm(vqT);
    eul      = WBM.utilities.rotm2eul(R, sequence);
end
