function qs = quatslerp(q1, q2, t)
    % This is a modified version of the function 'slerp' of Sagi Dalyot.
    % Source: <https://mathworks.com/matlabcentral/fileexchange/11827-slerp/content/slerp.m>
    %
    % Further details about the calculation can be found in:
    %   [1] Quaternions, Interpolation and Animation, Erik B. Dam & M. Koch & M. Lillholm, Technical Report DIKU-TR-98/5, 1998,
    %       Department of Computer Science, University of Copenhagen, <http://web.mit.edu/2.998/www/QuaternionReport1.pdf>,
    %       pp. 42-48, eq. (6.13).
    %   [2] Geometric Tools Engine, Documentation: <http://www.geometrictools.com/Documentation/Quaternions.pdf>, p. 6-7, eq. (21).
    %   [3] Animating Rotation with Quaternion Curves, Ken Shoemake, International Conference on Computer Graphics and Interactive Techniques,
    %       Volume 19, Number 3, 1985, <http://run.usc.edu/cs520-s15/assign2/p245-shoemake.pdf>.
    if ( (size(q1,1) ~= 4) || (size(q2,1) ~= 4) )
        error('quatslerp: %s', WBM.wbmErrorMsg.WRONG_VEC_DIM);
    end
    if ( (t < 0) || (t > 1) )
        error('quatslerp: Value t is out of range [0,1]!');
    end
    epsilon = 1e-12; % min. value to treat a number as zero ...

    switch t
        case 0 % special cases:
            qs = q1;
            return
        case 1
            qs = q2;
            return
        otherwise
            % compute the cosine of the angle between the two quaternion vectors ...
            cs_ang = q1.'*q2;

            theta_0 = acos(cs_ang); % angle between input vectors q1 and q2
            theta   = theta_0*t;    % angle between q1 and result qs

            if ((1 - cs_ang) <= epsilon)
                % If the angle theta is very close to 0 degrees, then
                % calculate by linear interpolation to avoid divisions
                % close to 0.
                qs = q1*(1 - t) + q2*t; % = q1 + t*(q2 - q1)

            elseif ((1 + cs_ang) <= epsilon)
                % If the dot-product (cosine of the angle) is negative, i.e.
                % the angle theta is very close to 180 degrees the result is
                % undefined. Thus, the quaternions have opposite handed-ness
                % and Slerp does not take the shortest direction to rotate.
                % Fix by inversing (rotation by 90 degrees) one of the unit
                % quaternions.
                q1_90 = zeros(4,1);
                q1_90(1,1) =  q1(4,1);
                q1_90(2,1) = -q1(3,1);
                q1_90(3,1) =  q1(2,1);
                q1_90(4,1) = -q1(1,1);

                pi_hf = pi*0.5;
                qs    = q1*sin((1 - t)*pi_hf) + q1_90*sin(t*pi_hf);
            else
                qs = (q1*sin((1 - t)*theta) + q2*sin(t*theta)) / sin(theta);
            end
    end
end
