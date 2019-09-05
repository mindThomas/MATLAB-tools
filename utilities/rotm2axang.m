function axang = rotm2axang(rotm)
    if ( (size(rotm,1) ~= 3) || (size(rotm,2) ~= 3) )
        error('rotm2axang: %s', WBM.wbmErrorMsg.WRONG_MAT_DIM);
    end
    axang   = zeros(4,1);
    epsilon = 1e-12; % min. value to treat a number as zero ...

    %% Translate a given rotation matrix R into the corresponding axis-angle representation (u, theta).
    % For further details about the computation, see:
    %   [1] Technical Concepts: Orientation, Rotation, Velocity and Acceleration and the SRM, P. Berner, Version 2.0, 2008,
    %       <http://sedris.org/wg8home/Documents/WG80485.pdf>, pp. 32-33.
    %   [2] A Mathematical Introduction to Robotic Manipulation, Murray & Li & Sastry, CRC Press, 1994, p. 30, eq. (2.17) & (2.18).
    %   [3] Modelling and Control of Robot Manipulators, L. Sciavicco & B. Siciliano, 2nd Edition, Springer, 2008,
    %       p. 35, formula (2.25).
    %   [4] Introduction to Robotics: Mechanics and Control, John J. Craig, 3rd Edition, Pearson/Prentice Hall, 2005,
    %       pp. 47-48, eq. (2.81) & (2.82).
    tr = rotm(1,1) + rotm(2,2) + rotm(3,3);
    if (abs(tr - 3) <= epsilon) % tr = 3 --> theta = 0:
        % Null rotation --> singularity: The rotation matrix R is the identity matrix and the axis of rotation u is undefined.
        % By convention, set u to the default value (0, 0, 1) according to the ISO/IEC IS 19775-1:2013 standard of the Web3D Consortium.
        % See: <http://www.web3d.org/documents/specifications/19775-1/V3.3/Part01/fieldsDef.html#SFRotationAndMFRotation>
        axang(3,1) = 1;
    elseif (abs(tr + 1) <= epsilon) % tr = -1 --> theta = pi:
        if ( (rotm(1,1) > rotm(2,2)) && (rotm(1,1) > rotm(3,3)) )
            u = vertcat(rotm(1,1)+1, rotm(1,2), rotm(1,3));
        elseif (rotm(2,2) > rotm(3,3))
            u = vertcat(rotm(2,1), rotm(2,2)+1, rotm(2,3));
        else
            u = vertcat(rotm(3,1), rotm(3,2), rotm(3,3)+1);
        end
        n = u.'*u;
        axang(1:3,1) = u./sqrt(n); % normalize
        axang(4,1)   = pi;
    else % general case, tr ~= 3 and tr ~= -1:
        axang(4,1) = acos((tr - 1)*0.5); % rotation angle theta within the range (0, pi).
        n_inv = 1/(2*sin(axang(4,1)));
        % unit vector u:
        axang(1,1) = (rotm(3,2) - rotm(2,3))*n_inv;
        axang(2,1) = (rotm(1,3) - rotm(3,1))*n_inv;
        axang(3,1) = (rotm(2,1) - rotm(1,2))*n_inv;
    end
end
