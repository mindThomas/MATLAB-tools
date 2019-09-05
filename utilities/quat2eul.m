function eul = quat2eul(quat, sequence)
    if (size(quat,1) ~= 4)
        error('quat2eul: %s', WBM.wbmErrorMsg.WRONG_VEC_DIM);
    end

    if ~exist('sequence', 'var')
        % use the default axis sequence ...
        sequence = 'ZYX';
    end
    eul = zeros(3,1);

    qnorm = quat.'*quat;
    if (qnorm > 1)
        quat = quat./sqrt(qnorm); % normalize
    end
    % scalar part:
    q_0 = quat(1,1);
    % vector part:
    q_1 = quat(2,1);
    q_2 = quat(3,1);
    q_3 = quat(4,1);

    %% Convert the given quaternion to the corresponding Euler angles in dependency of the
    %  order of the specified axis rotation sequence:
    % For further details see:
    %   [1] MATLAB TOOLBOX FOR RIGID BODY KINEMATICS, Hanspeter Schaub & John L. Junkins,
    %       9th AAS/AIAA Astrodynamics Specialist Conference, AAS 99-139, 1999, <http://hanspeterschaub.info/Papers/mtb1.pdf>,
    %       pp. 4-6, equations (15)-(17) & (18)-(20) by applying the homogeneous equation (6) with table 2.1 on p. 24 in [3].
    %   [2] GitHub: ShoolCode/ASEN 5010-Spacecraft Attitude Dynamics and Control/AIAA Software (2nd)/Matlab Toolbox,
    %       <https://github.com/ZachDischner/SchoolCode/tree/master/ASEN 5010-Spacecraft Attitude Dynamics and Control/AIAA Software (2nd)/Matlab Toolbox/>
    %   [3] Optimal Spacecraft Rotational Maneuvers, John L. Junkins & James D. Turner, Elsevier, 1986, p. 28 - eq. (2.55), p. 24 - table 2.1.
    %   [4] Shuttle Program. Euler Angles, Quaternions, and Transformation Matrices Working Relationships, D. M. Henderson, NASA, Mission Planning and Analysis Division,
    %       N77-31234/6, 1977, <http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770024290.pdf>, p. 7 - eq. (15), pp. 25-27 - Appendix A - eq. (10) & (12).
    switch sequence
        case 'ZYX'
            %                          r12                           r11
            eul(1,1) = atan2(2*(q_1*q_2 + q_0*q_3), q_0*q_0 + q_1*q_1 - q_2*q_2 - q_3*q_3); % theta_z
            %                         -r13
            eul(2,1) = asin(-2*(q_1*q_3 - q_0*q_2));                                        % theta_y
            %                          r23                           r33
            eul(3,1) = atan2(2*(q_2*q_3 + q_0*q_1), q_0*q_0 - q_1*q_1 - q_2*q_2 + q_3*q_3); % theta_x
        case 'ZYZ'
            % alternative computation by simplification, see [2] and [1, p. 5]:
            t1 = atan2(q_1, q_2);
            t2 = atan2(q_3, q_0);

            eul(1,1) = t2 - t1;                         % theta_z1
            %        = atan2(r32/r31)
            eul(2,1) = 2*acos(sqrt(q_0*q_0 + q_3*q_3)); % theta_y
            %        = acos(r33)
            eul(3,1) = t2 + t1;                         % theta_z2
            %        = atan2(r23/-r13)
        otherwise
            error('quat2eul: %s', WBM.wbmErrorMsg.UNKNOWN_AXIS_SEQ);
    end
end
