function qnorm = quatnorm(quat)
    if (size(quat,1) ~= 4)
        error('quatnorm: %s', WBM.wbmErrorMsg.WRONG_VEC_DIM);
    end

    %% Compute the norm of a quaternion q:
    % Further details can be found in:
    %   [1] GPGPU Programming for Games and Science, David H. Eberly, CRC Press, 2014, p. 296.
    %   [2] Geometric Tools Engine, Documentation: <http://www.geometrictools.com/Documentation/Quaternions.pdf>, p. 2, eq. (4).
    %   [3] Theory of Applied Robotics: Kinematics, Dynamics, and Control, Reza N. Jazar, 2nd Edition, Springer, 2010, p. 113, eq. (3.171).
    qnorm = quat.'*quat;
end
