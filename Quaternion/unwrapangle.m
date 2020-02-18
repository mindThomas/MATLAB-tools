function out = unwrapangle(angle, varargin)
    threshold = pi;
    
    if (nargin == 2)
        modulo = varargin{1};
    else
        modulo = 2*pi;
    end
    
    delta = diff(angle);
    jumps = find(abs(delta) > 0.75*modulo);
      
    offsets = [0; -modulo * cumsum(1-2*(delta(jumps)<0))];
    jumps = [1; jumps; length(angle)-1];  
    
    out = angle;
    for (i = 1:(length(jumps)-1))
        out(jumps(i)+1:jumps(i+1)) = out(jumps(i)+1:jumps(i+1)) + offsets(i);
    end
    
end
