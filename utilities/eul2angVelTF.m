function B = eul2angVelTF(eul, sequence)
    if (size(eul,1) ~= 3)
        error('eul2angVelTF: %s', WBM.wbmErrorMsg.WRONG_VEC_DIM);
    end

    if ~exist('sequence', 'var')
        % use the default sequence ...
        sequence = 'ZYX';
    end
    B = zeros(3,3);

    s_2 = sin(eul(2,1)); % sin(theta_y)
    c_2 = cos(eul(2,1)); % cos(theta_y)
    s_3 = sin(eul(3,1)); % sin(theta_x) or sin(theta_z2)
    c_3 = cos(eul(3,1)); % cos(theta_x) or cos(theta_z2)

    %% Body angle velocity transformations:
    %  They relating the time derivative of the Euler angle vector Theta to
    %  the body angular velocity vector w, s.t.
    %
    %       w = B(Theta)*dTheta/dt.
    %
    % Sources:
    %   [1] Optimal Spacecraft Rotational Maneuvers, John L. Junkins & James D. Turner, Elsevier, 1986, p. 21-24.
    %   [2] MATLAB TOOLBOX FOR RIGID BODY KINEMATICS, Hanspeter Schaub & John L. Junkins,
    %       9th AAS/AIAA Astrodynamics Specialist Conference, AAS 99-139, 1999, <http://hanspeterschaub.info/Papers/mtb1.pdf>, p. 12.
    %   [3] GitHub: ShoolCode/ASEN 5010-Spacecraft Attitude Dynamics and Control/AIAA Software (2nd)/Matlab Toolbox,
    %       <https://github.com/ZachDischner/SchoolCode/tree/master/ASEN 5010-Spacecraft Attitude Dynamics and Control/AIAA Software (2nd)/Matlab Toolbox/>
    %   [4] Rigid Body Kinematics and C++ Code, Sergio Pissanetzky, Scientific Controls, 2005, pp. 60-63.
    switch sequence
        case 'ZYX'
            %                |-s_2        0     1|
            % B(Theta_321) = | c_2*s_3    c_3   0|
            %                | c_2*c_3   -s_3   0|
            B(1,1) = -s_2;
            B(1,3) =  1;
            B(2,1) =  c_2*s_3;
            B(2,2) =  c_3;
            B(3,1) =  c_2*c_3;
            B(3,2) = -s_3;
        case 'ZYZ'
            %                |-s_2*c_3   s_3   0|
            % B(Theta_323) = | s_2*s_3   c_3   0|
            %                | c_2       0     1|
            B(1,1) = -s_2*c_3;
            B(1,2) =  s_3;
            B(2,1) =  s_2*s_3;
            B(2,2) =  c_3;
            B(3,1) =  c_2;
            B(3,3) =  1;
        otherwise
            error('eul2angVelTF: %s', WBM.wbmErrorMsg.UNKNOWN_AXIS_SEQ);
    end
end
