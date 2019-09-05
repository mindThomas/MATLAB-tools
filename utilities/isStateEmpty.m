function result = isStateEmpty(stParams)
    if ( ~isa(stParams, 'WBM.wbmStateParams') && ~isstruct(stParams) )
        error('isStateEmpty: %s', WBM.wbmErrorMsg.WRONG_DATA_TYPE);
    end
    result = false;

    if ( isempty(stParams.x_b) && isempty(stParams.qt_b) && isempty(stParams.q_j) && ...
         isempty(stParams.dx_b) && isempty(stParams.omega_b) && isempty(stParams.dq_j) )
        result = true;
    end
end
