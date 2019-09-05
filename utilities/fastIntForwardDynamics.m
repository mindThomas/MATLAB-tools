function [t, stmChi] = fastIntForwardDynamics(fhTrqControl, tspan, stvChi_0, robot_model, robot_config, ode_opt)
    if ~isa(fhTrqControl, 'function_handle')
        error('fastIntForwardDynamics: %s', WBM.wbmErrorMsg.WRONG_DATA_TYPE)
    end
    if ( (robot_config.stvLen == 0) || (robot_config.nCstrs == 0) )
        error('fastIntForwardDynamics: %s', WBM.wbmErrorMsg.VALUE_LTE_ZERO);
    end
    if ( isempty(robot_model.frict_coeff.v) || isempty(robot_model.frict_coeff.c) )
        error('fastIntForwardDynamics: %s', WBM.wbmErrorMsg.EMPTY_VECTOR);
    end

    if ~exist('ode_opt', 'var')
        % setup the default error tolerances ...
        ode_opt = odeset('RelTol', 1e-2, 'AbsTol', 1e-4);
    end

    fhFwdDyn    = @(t, chi)WBM.utilities.fastForwardDynamics(t, chi, fhTrqControl, robot_model, robot_config);
    [t, stmChi] = ode15s(fhFwdDyn, tspan, stvChi_0, ode_opt); % ODE-Solver
end
