function stvChi = params2stateVec(stParams)
    if WBM.utilities.isStateEmpty(stParams)
        error('params2stateVec: %s', WBM.wbmErrorMsg.EMPTY_DATA_TYPE);
    end
    if ~iscolumn(stParams.x_b)
        error('params2stateVec: %s', WBM.wbmErrorMsg.WRONG_DATA_TYPE);
    end

    stvChi = vertcat(stParams.x_b, stParams.qt_b, stParams.q_j, ...
                     stParams.dx_b, stParams.omega_b, stParams.dq_j);
end
