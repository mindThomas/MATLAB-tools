function dcm = quat2rotm(quat)
    if (length(quat) ~= 4)
        error('quat2rotm: %s', WBM.wbmErrorMsg.WRONG_VEC_DIM);
    end
    dcm = zeros(3,3);

    qnorm = norm(quat);
    if (qnorm > 1)
        quat = quat./qnorm; % normalize
    end
    % scalar (real) part:
    q_0 = quat(1);
    % vector (imaginary) part:
    q_1 = quat(2);
    q_2 = quat(3);
    q_3 = quat(4);

    %% Compute the Direction Cosine Matrix (DCM) by applying the Euler-Rodrigues Parameterization
    % for efficent computing of the rotation matrix R(s,r) in SO(3) with,
    %
    %       R(s,r) = I + (2*s)*S_r + 2*(S_r*S_r),
    %
    % where s denotes the scalar part of the given quaternion, S_r is the skew-symmetric matrix of
    % the vector part r of the quaternion, and I is a (3x3) identity matrix.
    % For more details about the parametrization formula, please check:
    %   [1] CONTRIBUTIONS TO THE AUTOMATIC CONTROL OF AERIAL VEHICLES, Minh Duc HUA, PhD-Thesis, 2009,
    %       <https://www-sop.inria.fr/act_recherche/formulaire/uploads/phd-425.pdf>, p. 101, eq. (3.7) & (3.8).
    %   [2] Lecture notes for 'Multibody System Dynamics', Prof. Bottasso, Aerospace Engineering Departement, Politecnico di Milano,
    %       <https://home.aero.polimi.it/trainelli/downloads/Bottasso_ThreeDimensionalRotations.pdf>, p. 33, eq. (1.160).
    %   [3] The Vectorial Parameterization of Rotation, Olivier A. Bauchau & Lorenzo Trainelli, Nonlinear Dynamics, 2003,
    %       <http://soliton.ae.gatech.edu/people/obauchau/publications/Bauchau+Trainelli03.pdf>, p. 16, eq. (51).
    %   [4] Theory of Applied Robotics: Kinematics, Dynamics, and Control, Reza N. Jazar, 2nd Edition, Springer, 2010, p. 103, eq. (3.82).
    %
    % Note: The direct assignment is computationally faster than using directly the formula above. Furthermore, if the above formula
    %       will be used, the matrix multiplications are causing accumulative round-off errors in the main diagonal of the DCM-matrix
    %       (around ??0.7*1e-03). To be more accurate, it is better to use the derived transformation matrix of the formula.
    dcm(1,1) = q_0*q_0 + q_1*q_1 - q_2*q_2 - q_3*q_3;
    dcm(1,2) = 2*(q_1*q_2 - q_0*q_3);
    dcm(1,3) = 2*(q_0*q_2 + q_1*q_3);

    dcm(2,1) = 2*(q_0*q_3 + q_1*q_2);
    dcm(2,2) = q_0*q_0 - q_1*q_1 + q_2*q_2 - q_3*q_3;
    dcm(2,3) = 2*(q_2*q_3 - q_0*q_1);

    dcm(3,1) = 2*(q_1*q_3 - q_0*q_2);
    dcm(3,2) = 2*(q_0*q_1 + q_2*q_3);
    dcm(3,3) = q_0*q_0 - q_1*q_1 - q_2*q_2 + q_3*q_3;
end
