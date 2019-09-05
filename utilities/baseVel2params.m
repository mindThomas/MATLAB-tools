function [dx_b, omega_b] = baseVel2params(v_b)
    if iscolumn(v_b)
        if (size(v_b,1) ~= 6)
           error('baseVel2params: %s', WBM.wbmErrorMsg.WRONG_VEC_DIM);
        end

        dx_b    = v_b(1:3,1);
        omega_b = v_b(4:6,1);
        return
    elseif ismatrix(v_b)
        [nRows, nCols] = size(v_b);
        if (nCols ~= 6)
            error('baseVel2params: %s', WBM.wbmErrorMsg.WRONG_MAT_DIM);
        end

        dx_b    = v_b(1:nRows,1:3);
        omega_b = v_b(1:nRows,4:6);
        return
    end
    % else ...
    error('baseVel2params: %s', WBM.wbmErrorMsg.WRONG_DATA_TYPE);
end
