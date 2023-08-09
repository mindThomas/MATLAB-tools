%
% principalValue(angle, min) maps angle into interval [min, min + 2*PI)
% the implementation is only intended to be used for "small" values of angle, i.e.
% |angle| < 1e2
%
% typically, min == -PI, or min == 0.
%
% if min is not specified the value min == -PI is taken
% and this function thus maps between [-PI, PI)
%
function angle = principalValue(angle, varargin)   
    if (length(angle) > 1)
        angle = angle(1);
    end
    
    if (length(varargin) == 2)
        minimum = varargin{1};
        tau = varargin{2};
    elseif (length(varargin) == 1)
        minimum = varargin{1};
        tau = 2*pi;
    else
        minimum = -pi;
        tau = 2*pi;
    end

    while (angle < minimum)
        angle = angle + tau;
    end

    while (minimum + tau <= angle)
        angle = angle - tau;
    end
end