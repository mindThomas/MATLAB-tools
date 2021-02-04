function heading = rotm2heading(rotm, varargin)
    eul = rotm2eul(rotm, 'ZYX');
    heading = eul(1);
end
