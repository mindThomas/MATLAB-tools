function quat = quatmult(q1, q2)
    if ( (size(q1,1) ~= 4) || (size(q2,1) ~= 4) )
        error('quatmult: %s', WBM.wbmErrorMsg.WRONG_VEC_DIM);
    end
    quat = zeros(4,1);

    %% Compute the product of two quaternions q1 and q2 with q = q1*q2:
    %  Note: The quaternion multiplication is in general not commutative.
    %   scalar part:
    %   q_0  =  q1_0*q2_0  -  q1_1*q2_1  -  q1_2*q2_2  -  q1_3*q2_3
    %   vector part:
    %   q_1  =  q1_0*q2_1  +  q2_0*q1_1  +  q1_2*q2_3  -  q1_3*q2_2
    %   q_2  =  q1_0*q2_2  +  q2_0*q1_2  +  q1_3*q2_1  -  q1_1*q2_3
    %   q_3  =  q1_0*q2_3  +  q2_0*q1_3  +  q1_1*q2_2  -  q1_2*q2_1
    %
    % For further details about the quaternion multiplcation, see:
    %   [1] GPGPU Programming for Games and Science, David H. Eberly, CRC Press, 2014, p. 295, eq. (6.78).
    %   [2] Theory of Applied Robotics: Kinematics, Dynamics, and Control, Reza N. Jazar, 2nd Edition, Springer, 2010, p. 112, eq. (3.166).
    %   [3] Optimal Spacecraft Rotational Maneuvers, John L. Junkins & James D. Turner, Elsevier, 1986, p. 38, eq. (2.84).
    %   [4] Geometric Tools Engine, Documentation: <http://www.geometrictools.com/Documentation/Quaternions.pdf>, p. 2, eq. (2).
    %   [5] 3D Game Engine Programming, Understanding Quaternions: <http://www.3dgep.com/understanding-quaternions/#Quaternion_Products>.
    % scalar part:
    quat(1,1) = q1(1,1)*q2(1,1) - q1(2,1)*q2(2,1) - q1(3,1)*q2(3,1) - q1(4,1)*q2(4,1);
    % vector part:
    quat(2,1) = q1(1,1)*q2(2,1) + q2(1,1)*q1(2,1) + q1(3,1)*q2(4,1) - q1(4,1)*q2(3,1);
    quat(3,1) = q1(1,1)*q2(3,1) + q2(1,1)*q1(3,1) + q1(4,1)*q2(2,1) - q1(2,1)*q2(4,1);
    quat(4,1) = q1(1,1)*q2(4,1) + q2(1,1)*q1(4,1) + q1(2,1)*q2(3,1) - q1(3,1)*q2(2,1);
end
