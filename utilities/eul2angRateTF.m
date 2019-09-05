function B_inv = eul2angRateTF(eul, sequence)
    if (size(eul,1) ~= 3)
        error('eul2angRateTF: %s', WBM.wbmErrorMsg.WRONG_VEC_DIM);
    end

    if ~exist('sequence', 'var')
        % use the default sequence ...
        sequence = 'ZYX';
    end
    B_inv = zeros(3,3);

    s_2 = sin(eul(2,1)); % sin(theta_y)
    c_2 = cos(eul(2,1)); % cos(theta_y)
    s_3 = sin(eul(3,1)); % sin(theta_x) or sin(theta_z2)
    c_3 = cos(eul(3,1)); % cos(theta_x) or cos(theta_z2)

    %% Euler angle rate transformation matrices:
    %  They relating the angular velocity vector w of the body to the
    %  time derivative of the Euler angle vector Theta, s.t.
    %
    %       dTheta/dt = B(Theta)^(-1)*w.
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
            %                             |0     s_3           c_3|
            % B(Theta_321)^(-1) = 1/c_2 * |0     c_2*c_3  -c_2*s_3|
            %                             |c_2   s_2*s_3   s_2*c_3|
            B_inv(1,2) = s_3;
            B_inv(1,3) = c_3;
            B_inv(2,2) = c_2*c_3;
            B_inv(2,3) = -c_2*s_3;
            B_inv(3,1) = c_2;
            B_inv(3,2) = s_2*s_3;
            B_inv(3,3) = s_2*c_3;
            B_inv = B_inv/c_2;
        case 'ZYZ'
            %                             |-c_3        s_3         0|
            % B(Theta_323)^(-1) = 1/s_2 * | s_2*s_3    s_2*c_3     0|
            %                             | c_2*c_3   -c_2*s_3   s_2|
            B_inv(1,1) = -c_3;
            B_inv(1,2) = s_3;
            B_inv(2,1) = s_2*s_3;
            B_inv(2,2) = s_2*c_3;
            B_inv(3,1) = c_2*c_3;
            B_inv(3,2) = -c_2*s_3;
            B_inv(3,3) = s_2;
            B_inv = B_inv/s_2;
        otherwise
            error('eul2angRateTF: %s', WBM.wbmErrorMsg.UNKNOWN_AXIS_SEQ);
    end
end
