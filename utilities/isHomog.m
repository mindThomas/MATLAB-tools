function result = isHomog(tform, epsilon)
    if ( (size(tform,1) ~= 4) || (size(tform,2) ~= 4) )
        error('isHomog: %s', WBM.wbmErrorMsg.WRONG_MAT_DIM);
    end
    result = false;

    if ~exist('epsilon', 'var')
        epsilon = 1e-12; % min. value to treat a number as zero ...
    end

    if (abs(det(tform) - 1) <= epsilon)
        result = true;
    end
end
