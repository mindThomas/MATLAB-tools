function tf = rotm2tf(rotm, varargin)
    origin = zeros(3,1);
    if (length(varargin) == 1 && length(varargin{1}) == 3)
        origin = varargin{1};
    end
    tf = [rotm, origin; 0,0,0,1];
end
