function quat = eul2quat(eul, sequence)
    if (size(eul,1) ~= 3)
        error('eul2quat: %s', WBM.wbmErrorMsg.WRONG_VEC_DIM);
    end

    if ~exist('sequence', 'var')
        % use the default axis sequence ...
        sequence = 'ZYX';
    end
    quat = zeros(4,1);

    %% Translate the given Euler angle with a specified axis rotation sequence into the corresponding quaternion:
    % For further details see:
    %   [1] MATLAB TOOLBOX FOR RIGID BODY KINEMATICS, Hanspeter Schaub & John L. Junkins,
    %       9th AAS/AIAA Astrodynamics Specialist Conference, AAS 99-139, 1999, <http://hanspeterschaub.info/Papers/mtb1.pdf>,
    %       p. 7, equations (26a)-(26d).
    %   [2] GitHub: ShoolCode/ASEN 5010-Spacecraft Attitude Dynamics and Control/AIAA Software (2nd)/Matlab Toolbox,
    %       <https://github.com/ZachDischner/SchoolCode/tree/master/ASEN 5010-Spacecraft Attitude Dynamics and Control/AIAA Software (2nd)/Matlab Toolbox/>
    %   [3] Optimal Spacecraft Rotational Maneuvers, John L. Junkins & James D. Turner, Elsevier, 1986, pp. 31-34, table 2.2.
    %   [4] Shuttle Program. Euler Angles, Quaternions, and Transformation Matrices Working Relationships, D. M. Henderson, NASA, Mission Planning and Analysis Division,
    %       N77-31234/6, 1977, <http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770024290.pdf>, pp. 25-27, Appendix A - eq. (10) & (12).
    switch sequence
        case 'ZYX'
			c_1 = cos(eul(1,1)*0.5);
			s_1 = sin(eul(1,1)*0.5);
			c_2 = cos(eul(2,1)*0.5);
			s_2 = sin(eul(2,1)*0.5);
			c_3 = cos(eul(3,1)*0.5);
			s_3 = sin(eul(3,1)*0.5);

			quat(1,1) = c_1*c_2*c_3 + s_1*s_2*s_3;
			quat(2,1) = c_1*c_2*s_3 - s_1*s_2*c_3;
			quat(3,1) = c_1*s_2*c_3 + s_1*c_2*s_3;
			quat(4,1) = s_1*c_2*c_3 - c_1*s_2*s_3;
        case 'ZYZ'
			t_1 = eul(1,1)*0.5;
			t_2 = eul(2,1)*0.5;
			t_3 = eul(3,1)*0.5;

			quat(1,1) = cos(t_2)*cos( t_1 + t_3);
			quat(2,1) = sin(t_2)*sin(-t_1 + t_3);
			quat(3,1) = sin(t_2)*cos(-t_1 + t_3);
			quat(4,1) = cos(t_2)*sin( t_1 + t_3);
        otherwise
            error('eul2quat: %s', WBM.wbmErrorMsg.UNKNOWN_AXIS_SEQ);
    end
end
