function result = isRadians(x)
    % A modified version of the function 'isradians' of Michaël Zugaro which is under GPL version 3.
    % Source: <http://fmatoolbox.sourceforge.net/API/FMAToolbox/Helpers/isradians.html>
    result = 0;

    if ~isa(x, 'double')
        return
    end

    % ignore NaN ...
    x = x(~isnan(x));
    if isempty(x)
        return
    end

    if ~isscalar(x)
        % find the min. and the max. value of the array ...
        minv = min(x);
        maxv = max(x);
    else % x is a single value ...
        minv = x;
        maxv = minv;
    end

    % range test ...
    if ( (minv >= -pi) && (maxv <= pi) )
        result = 1;
        return
    end
    if ( (minv >= 0) && (maxv <= 2*pi) )
        result = 2;
    end
end
