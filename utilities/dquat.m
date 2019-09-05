function dq = dquat(quat, omega, varargin)
    if (size(quat,1) ~= 4)
        error('dquat: %s', WBM.wbmErrorMsg.WRONG_VEC_DIM);
    end
    if (size(omega,1) ~= 3)
        error('dquat: %s', WBM.wbmErrorMsg.WRONG_VEC_DIM);
    end

    switch nargin
        case 2
            % default values:
            k      = 1;     % gain (drives the norm of the quaternion to 1, ε should become nonzero)
            chktol = false; % tolerance check is deactivated
        case 3
            if islogical(varargin{1,1})
                k      = 1; % use default
                chktol = varargin{1,1};
            elseif isreal(varargin{1,1})
                k      = varargin{1,1};
                chktol = false; % use default
            else
                error('dquat: %s', WBM.wbmErrorMsg.WRONG_DATA_TYPE);
            end
        case 4
            k      = varargin{1,1};
            chktol = varargin{1,2};
        otherwise
            error('dquat: %s', WBM.wbmErrorMsg.WRONG_ARG);
    end

    qnorm = quat.'*quat;
    epsilon = 1 - sqrt(qnorm); % = 1 - norm(quat)

    if chktol
        eps_abs = abs(epsilon);
        % check tolerances:
        if (eps_abs >= 1e-4)
            if (eps_abs > 0.5)
                error('dquat: Quaternion left SO(3), abs(epsilon) > 0.5!');
            end
            % else ...
            warning('dquat: Quaternion is leaving SO(3), abs(epsilon) >= 1e-4.');
        end
    end
    % Create the conjugate product matrix "Omega" of the angular velocity omega:
    % Further details about the product matrix (operator) Omega can be taken from:
    %   [1] CONTRIBUTIONS TO THE AUTOMATIC CONTROL OF AERIAL VEHICLES, Minh Duc HUA, PhD-Thesis, 2009,
    %       <https://www-sop.inria.fr/act_recherche/formulaire/uploads/phd-425.pdf>, p. 101.
    %   [2] The Vectorial Parameterization of Rotation, Olivier A. Bauchau & Lorenzo Trainelli, Nonlinear Dynamics, 2003,
    %       <http://soliton.ae.gatech.edu/people/obauchau/publications/Bauchau+Trainelli03.pdf>, p. 16-18, formulas (51) & (73).
    %   [3] Quaternion kinematics for the error-state KF, Joan Solà, Universitat Politècnica de Catalunya, 2016,
    %       <http://www.iri.upc.edu/people/jsola/JoanSola/objectes/notes/kinematics.pdf>, p. 19, formulas (97) & (98).
    Omega_x = zeros(4,4);
    Omega_x(2:4,2:4) = -WBM.utilities.skewm(omega); % (skew-symmetric) cross product matrix
    Omega_x(1,2:4)   = -omega.';
    Omega_x(2:4,1)   =  omega;

    % This calculation of the quaternion derivative has additionally a
    % computational "hack" (k*ε*q) to compute the derivative in the case
    % that the quaternion is not normalized, s.t. the vector sum is still
    % equal to 1. The factor ε will be 0 if the quaternion is unitary.
    % For further details see:
    %   <http://mathworks.com/help/aeroblks/customvariablemass6dofquaternion.html>
    %   and <https://www.physicsforums.com/threads/quaternion-derivative.706475>
    dq = 0.5*(Omega_x*quat) + (k*epsilon)*quat;
end
