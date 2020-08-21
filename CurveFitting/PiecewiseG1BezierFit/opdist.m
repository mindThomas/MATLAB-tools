function se2 = opdist(id,Q,P,ang)

% function se2 = opdist(id,Q,P,ang). This function is the obj
% ective function to be minimized during segment optimization.
% It receives subcomponents from a curve's composite vector, a
% segment at a time. It returns the sum error for a segment.
% It was written by E. J. Lane.

n = length(Q);
C = ctpts(P,ang,id); 	% Call to compute the segments
% control points.
se2 = 0;

% Loop to find distance error in a segment.
for j = 1:n
    np = NearestPoint(C', Q(:,j)');
    
    if (j==1) && isequal(np, C(:,1)')
        d = zeros(1,2);
    elseif (j==n) && isequal(np, C(:,4)')
        d = zeros(1,2);
    else
        d = (Q(:,j)' - np);
    end
    se2 = se2 + d*d';
end