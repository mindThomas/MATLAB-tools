function stParams = fastGetStateParams(stvChi, stvLen, ndof)
    if ~iscolumn(stvChi)
       error('fastGetStateParams: %s', WBM.wbmErrorMsg.WRONG_VEC_DIM);
    end

    % get the base/joint positions and the base orientation ...
    stParams.x_b  = stvChi(1:3,1);
    stParams.qt_b = stvChi(4:7,1);
    stParams.q_j  = stvChi(8:ndof+7,1);
    % get the velocities ...
    stParams.dx_b    = stvChi(ndof+8:ndof+10,1);
    stParams.omega_b = stvChi(ndof+11:ndof+13,1);
    stParams.dq_j    = stvChi(ndof+14:stvLen,1);
end
