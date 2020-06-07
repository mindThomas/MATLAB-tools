function T = trotz(angle, varargin)
    if (length(varargin) == 1)
        t = varargin{1};
    else
        t = [0;0;0];
    end
    
    T = [rotz(angle), t; 0,0,0,1];
end