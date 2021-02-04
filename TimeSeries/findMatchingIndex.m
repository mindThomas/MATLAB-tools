% This assumes a sorted time index list
function idx = findMatchingIndex(time, tSearch, varargin)
    if (nargin == 3)
        idx = varargin{1};
    else
        idx = 1;
    end
    
    n = length(time);
    if (idx < 1)
        idx = 1;
    end
    if (idx > n)
        idx = n;
    end
        
    while (idx < (n-1))
        if (time(idx) <= tSearch && time(idx+1) > tSearch)     
            if (abs(time(idx) - tSearch) > abs(time(idx+1) - tSearch))
                idx = idx + 1;
            end
            return;
        end
        idx = idx + 1;
    end
    
    % Error, time was not found
    idx = -1;    
end