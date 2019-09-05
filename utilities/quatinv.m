function q_inv = quatinv(quat)
    q_conj  = WBM.utilities.quatconj(quat);
    epsilon = 1e-12; % min. value to treat a number as zero ...

    %% Compute the inverse q^(-1) of a quaternion q:
    % Further details can be found in:
    %   [1] GPGPU Programming for Games and Science, David H. Eberly, CRC Press, 2014, p. 296.
    %   [2] Geometric Tools Engine, Documentation: <http://www.geometrictools.com/Documentation/Quaternions.pdf>, p. 2, eq. (5).
    %   [3] Theory of Applied Robotics: Kinematics, Dynamics, and Control, Reza N. Jazar, 2nd Edition, Springer, 2010, p. 113, eq. (3.173).
    qnorm = quat.'*quat;
    if ((qnorm - 1) <= epsilon)
        q_inv = q_conj;
        return
    end
    % else ...
    qn_inv = 1/qnorm;
    q_inv  = q_conj*qn_inv;
end
