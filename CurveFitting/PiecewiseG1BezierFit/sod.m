function sumd = sod(C,Q,dpkpc)

% function sumofdist = sod(C,Q,dpkpc). This function receives inputs:
% control points,C, data points,Q, and the dividing points or knot
% sequence. It finds the closest point on the curve for a given data
% point and computes the distance error. The function returns the
% sum of the distance squared from the data points to their nearest
% point on the curve segment. This was written by M. R. Holmes.

n     = length(C);
[r,s] = size(Q);
y     = dpkpc;
cntr  = 0;
sum   = 0;
for i = 1:3:n-3
    cntr = cntr + 1;
    % Loop to find distances from
    % data points to closest point
    % on the curve.
    for j = y(cntr):y(cntr+1)
        np = NearestPoint(C(:,i:i+3)', Q(:,j)');
        d = (Q(:,j)'- np);
        sum = sum + d * d'; % Note: NearestPoint is a MATLAB
        if (j==y(cntr)) && (i>1) 	% interface program written by
            d2  = d*d'; 		% Dr C. Borges for some "C" rou-
            ds2 = ds*ds'; 		% tines obtained from 'Solving
            dm  = max(d2,ds2); 	% the Nearest Point-on-Curve
            sum = sum - dm;     % Problem' and 'A Bezier Curve
                                % Root-Finder' developed P.J.
                                % Schneider in "Graphics Gems"
                                % Academic Press, 1990.
        end
    end
    ds = d;
end
sumd = sum;