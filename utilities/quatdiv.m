function quat = quatdiv(q1, q2)
    if ( (size(q1,1) ~= 4) || (size(q2,1) ~= 4) )
        error('quatdiv: %s', WBM.wbmErrorMsg.WRONG_VEC_DIM);
    end
    quat   = zeros(4,1);
    qn_inv = 1/(q2.'*q2); % the .' is crucial to avoid computing the conjugate!

    %% Calculate the division of two quaternions q1 and q2 with q = q1/q2.
    % The division can be calculated by multiplication of the inverse.
    % Since the quaternion multiplication is not commutative, we have to
    % distinguish between 'premultiplikation' with q = q2^(-1)*q1, and
    % 'postmultiplikation' with q = q1*q2^(-1).
    % Note: This function uses the premultiplication q = q2^(-1)*q1.
    % This can be done either by calling,
    %
    %   q2_inv = quatinv(q2);
    %   quat   = quatmult(q2_inv, q1);
    %
    % or alternatively, directly in transformation form (which is faster):
    %   scalar part:
    %   q_0  =  (q2_0*q1_0  +  q2_1*q1_1  +  q2_2*q1_2  +  q2_3*q1_3) / (q2_0^2 + q2_1^2 + q2_2^2 + q2_3^2)
    %   vector part:
    %   q_1  =  (q2_0*q1_1  -  q2_1*q1_0  -  q2_2*q1_3  +  q2_3*q1_2) / (q2_0^2 + q2_1^2 + q2_2^2 + q2_3^2)
    %   q_2  =  (q2_0*q1_2  +  q2_1*q1_3  -  q2_2*q1_0  -  q2_3*q1_1) / (q2_0^2 + q2_1^2 + q2_2^2 + q2_3^2)
    %   q_3  =  (q2_0*q1_3  -  q2_1*q1_2  +  q2_2*q1_1  -  q2_3*q1_0) / (q2_0^2 + q2_1^2 + q2_2^2 + q2_3^2)
    %
    % Furter details about the quaternion division can be taken from:
    %   [1] MathWorks, Documentation - Aerospace Toolbox: <http://mathworks.com/help/aerotbx/ug/quatdivide.html>.
    %   [2] Geometric Tools Engine, Documentation: <http://www.geometrictools.com/Documentation/Quaternions.pdf>, p. 2.
    %   [3] Aircraft Control and Simulation: Dynamics, Controls Design, and Autonomous Systems, Stevens & Lewis & Johnson,
    %       3rd Edition, Wiley-Blackwell, 2015, pp. 46-48.
    % scalar part:
    quat(1,1) = (q2(1,1)*q1(1,1) + q2(2,1)*q1(2,1) + q2(3,1)*q1(3,1) + q2(4,1)*q1(4,1))*qn_inv;
    % vector part:
    quat(2,1) = (q2(1,1)*q1(2,1) - q2(2,1)*q1(1,1) - q2(3,1)*q1(4,1) + q2(4,1)*q1(3,1))*qn_inv;
    quat(3,1) = (q2(1,1)*q1(3,1) + q2(2,1)*q1(4,1) - q2(3,1)*q1(1,1) - q2(4,1)*q1(2,1))*qn_inv;
    quat(4,1) = (q2(1,1)*q1(4,1) - q2(2,1)*q1(3,1) + q2(3,1)*q1(2,1) - q2(4,1)*q1(1,1))*qn_inv;
end
