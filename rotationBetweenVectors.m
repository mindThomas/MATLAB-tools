function R = rotationBetweenVectors(v1, v2)
    
    % Compute the rotation that v1 into v2
    %
    % If the two vectors are direction vectors in two frame, v1 in frame A
    % and v2 in frame B, the computed rotation matrix will rotate points
    % from frame B into frame A   

%     % Use cross product to compute rotation vector
%     r = cross(v1, v2);
%     % Use dot product to compute rotation amount
%     theta = acos(dot(v1, v2) / (norm(v1)*norm(v2)))
%     
%     % Construct quaternion from axis-angle theorem
%     r_ = r / norm(r);    
%     q = [cos(theta/2); r_ * sin(theta/2)];
%     
%     % Convert quaternion to rotation matrix
%     R = quat2rotm(q)
    
    %% Using Rodrigues' rotation formula   
    % Normalize vectors
    v1_ = v1 / norm(v1);
    v2_ = v2 / norm(v2);
    
    r = cross(v1_, v2_);
    cos_theta = dot(v1_, v2_);
    sin_theta = norm(r);
    theta = asin(sin_theta)
    
    R = eye(3) + skew(r) + skew(r)*skew(r)*1/(1+cos_theta)
    
end