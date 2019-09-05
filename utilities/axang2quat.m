function quat = axang2quat(axang)
    if (size(axang,1) ~= 4)
        error('axang2quat: %s', WBM.wbmErrorMsg.WRONG_VEC_DIM);
    end
    quat = zeros(4,1);
    u    = axang(1:3,1); % rotation axis vector

    n = u.'*u; % the .' is crucial to avoid computing the conjugate!
    if (n > 1)
        u = u./sqrt(n); % normalize u
    end

    %% Translate the axis-angle representation (u, theta) into the corresponding quaternion:
    % For further details, see:
    %   [1] MATLAB TOOLBOX FOR RIGID BODY KINEMATICS, Hanspeter Schaub & John L. Junkins,
    %       9th AAS/AIAA Astrodynamics Specialist Conference, AAS 99-139, 1999, <http://hanspeterschaub.info/Papers/mtb1.pdf>, p. 3, eq. (2).
    %   [2] GitHub: ShoolCode/ASEN 5010-Spacecraft Attitude Dynamics and Control/AIAA Software (2nd)/Matlab Toolbox,
    %       <https://github.com/ZachDischner/SchoolCode/tree/master/ASEN 5010-Spacecraft Attitude Dynamics and Control/AIAA Software (2nd)/Matlab Toolbox/>
    %   [3] Geometric Tools Engine, Documentation: <http://www.geometrictools.com/Documentation/RotationIssues.pdf>, p. 8.
    %   [4] Representing Attitude: Euler Angles, Unit Quaternions, and Rotation Vectors, J. Diebel, Stanford University, 2006,
    %       <https://www.astro.rug.nl/software/kapteyn/_downloads/attitude.pdf>, p. 17, eq. (175).
    ang_hf = axang(4,1)*0.5; % rotation angle theta / 2 (in radians)
    s_a = sin(ang_hf);

    quat(1,1) = cos(ang_hf); % scalar part
    quat(2,1) = u(1,1)*s_a;  % vector part
    quat(3,1) = u(2,1)*s_a;
    quat(4,1) = u(3,1)*s_a;
end
