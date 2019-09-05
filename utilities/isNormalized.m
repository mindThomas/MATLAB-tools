function result = isNormalized(v, epsilon)
    if ~iscolumn(v)
        error('isNormalized: %s', WBM.wbmErrorMsg.WRONG_VEC_DIM);
    end
    if ~exist('epsilon', 'var')
        epsilon = 1e-12; % min. value to treat a number as zero ...
    end
    result = false;

    if ((v.'*v - 1) <= epsilon) % the .' is crucial to avoid computing the conjugate!
        % the vector is already normalized ...
        result = true;
    end
end
