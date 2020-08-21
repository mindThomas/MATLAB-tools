function C = ctpts(P, ang, dt)

% C = ctpts(P,ang,dt). This function takes knot points,P; angles
% of the tangent vectors, ang; distances between successive knot
% points, dt; as input. It then computes the positions for the
% control points. It was written by M. R. Holmes.

n = length(P);
T = [];   % ECR

for k= 2:n-1
    % Converts the interior angles into their x and y components.
    u = [cos(ang(k)); sin(ang(k))];
    
    %Assembles the vector knot points with their
    %adjacent interior control points.
    T = [T P(:,k)-u*dt(2,k-1) P(:,k) P(:,k)+u*dt(1,k)];
end

u1 = [cos(ang(1)) ; sin(ang(1))]; % Converts the first and last
un = [cos(ang(n)) ; sin(ang(n))]; % angles into their x and y components.

% Assembles the vector of all control points
C  = [P(:,1) P(:,1)+u1*dt(1,1) T P(:,n)-un*dt(2,n-1) P(:,n)];