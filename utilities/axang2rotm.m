function dcm = axang2rotm(axang)
    if (size(axang,1) ~= 4)
        error('axang2rotm: %s', WBM.wbmErrorMsg.WRONG_VEC_DIM);
    end
    dcm = zeros(3,3);
    u   = axang(1:3,1); % rotation axis vector

    n = u.'*u;
    if (n > 1)
        u = u./sqrt(n); % normalize u
    end
    % axis elements:
    u_1    = u(1,1);
    u_2    = u(2,1);
    u_3    = u(3,1);
    % angles:
    c_t    = cos(axang(4,1));
    s_t    = sin(axang(4,1));
    vers_t = 1 - c_t; % = versine(theta) = 2*sin^2(theta/2)

    %% Compute the Direction Cosine Matrix (DCM) by applying the Rodrigues' formula
    %
    %       R(u,theta) = I + sin(theta)*S_u + (1 - cos(theta))*(S_u*S_u),
    %
    % where S_u denotes the skew-symmetric matrix of the rotation axis u, and I = eye(3,3).
    % For further details about the Rodrigues' formula and the corresponding rotation matrix R(u,theta), see:
    %   [1] Theory of Applied Robotics: Kinematics, Dynamics, and Control, Reza N. Jazar, 2nd Edition, Springer, 2010, p. 92, eq. (3.4) & (3.5).
    %   [2] Rigid Body Kinematics and C++ Code, Sergio Pissanetzky, Scientific Controls, 2005, p. 14, eq. (1.47 & (1.48).
    %   [3] A Mathematical Introduction to Robotic Manipulation, Murray & Li & Sastry, CRC Press, 1994, p. 28, eq. (2.14).
    %   [4] Technical Concepts: Orientation, Rotation, Velocity and Acceleration and the SRM, P. Berner, Version 2.0, 2008,
    %       <http://sedris.org/wg8home/Documents/WG80485.pdf>, pp. 33-34, eq. (1.25) & (1.26).
    % Note: The direct assignment is computationally faster than applying directly the formula above.
    dcm(1,1) = u_1*u_1*vers_t + c_t;
    dcm(1,2) = u_1*u_2*vers_t - u_3*s_t;
    dcm(1,3) = u_1*u_3*vers_t + u_2*s_t;

    dcm(2,1) = u_1*u_2*vers_t + u_3*s_t;
    dcm(2,2) = u_2*u_2*vers_t + c_t;
    dcm(2,3) = u_2*u_3*vers_t - u_1*s_t;

    dcm(3,1) = u_1*u_3*vers_t - u_2*s_t;
    dcm(3,2) = u_2*u_3*vers_t + u_1*s_t;
    dcm(3,3) = u_3*u_3*vers_t + c_t;
end
