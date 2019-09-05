function dv_b = generalizedBaseAcc(M, c_qv, ddq_j, ndof)
    % GENERALIZEDBASEACC computes the generalized floating base acceleration for a hybrid-dynamic
    % system.
    %
    % Useful for example in inverse dynamics, when the base acceleration is unknown.
    %
    %   INPUT ARGUMENTS:
    %       M     -- ((nDoF+6) x (nDoF+6)) floating base mass matrix
    %       c_qv  -- ((nDoF+6) x 1) generalized bias force vector
    %       ddq_j -- (nDoF x 1) joint angle acceleration vector (rad/s^2)
    %       ndof  -- number of degrees of freedom of the robot (optional)
    %
    %   OUTPUT ARGUMENTS:
    %       dv_b -- (6 x 1) generalized base acceleration vector.
    %
    % Further details about the formula are available at:
    %   [1] Rigid Body Dynamics Algorithms, Roy Featherstone, Springer, 2008,
    %       chapter 9.3-9.5, pp. 180-184, eq. (9.13) & (9.24).
    %   [2] Informatics in Control, Automation and Robotics, J. A. Cetto & J. Ferrier & J. Filipe,
    %       Lecture Notes in Electrical Engineering, volume 89, Springer, 2011, p. 14, eq. (36) & (37).
    %   [3] Dynamics of Tree-Type Robotic Systems, S. V. Shah & S. K. Saha & J. K. Dutt,
    %       Intelligent Systems, Control and Automation: Science and Engineering, Volume 62, Springer, 2012,
    %       p. 119, eq. (7.1) & (7.2).
    %
    % Author: Martin Neururer (martin.neururer@gmail.com); Genova, Jan 2017
    if exist('ndof', 'var')
        n = ndof + 6;
    else
        n = size(M,2);
    end

    h_0  = c_qv(1:6,1);
    M_00 = M(1:6,1:6);
    M_01 = M(1:6,7:n);

    dv_b = -M_00 \ (M_01*ddq_j + h_0);
end
