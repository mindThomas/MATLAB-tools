function out = unwrapstep(angle, prev_angle, prev_out_angle, varargin)
    if (nargin == 4)
        modulo = varargin{1};
    else
        modulo = 2*pi;
    end
 
    delta = angle - prev_angle;
    if (abs(delta) > 0.8*modulo)
        delta = delta - modulo*floor((angle - prev_angle) / modulo);
    end
    out = prev_out_angle + delta;
end
