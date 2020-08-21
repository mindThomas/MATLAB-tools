function [D,lambda,isConvex] = findClosestPoint(A,B,C)
% [D,lambda,isConvex] = findClosestPoint(A,B,C)
%
% Given a line segment AB, and an arbitrary point C, compute the point D on
% the line AB which is closest to point C.
%
% INPUTS:
%   A = [2,n] = [Ax;Ay] = First point on the line
%   B = [2,n] = [Bx;By] = Second point on the line
%   C = [2,n] = [Cx;Cy] = Arbitrary point of interest
%
% OUTPUTS:
%   D = [2,n] = [Dx;Dy] = point on AB that is nearest to C
%   lambda = [1,n] = lagrange multiplier:   D = lambda*A + (1-lambda)*B
%   isConvex = true if D is a convex combination of A and B
%               isConvex =   lambda >=0  &&  lambda <= 1
%

Ax = A(1,:);
Ay = A(2,:);

Bx = B(1,:);
By = B(2,:);

Cx = C(1,:);
Cy = C(2,:);

CBx = Cx-Bx;
ABx = Ax-Bx;
CBy = Cy-By;
ABy = Ay-By;

% This is derived from the constraint that: 0 = dot(AB,CD)
lambda = (CBx.*ABx + CBy.*ABy)./( ABx.*ABx + ABy.*ABy );
isConvex = lambda >= 0  & lambda <= 1;

Dx = lambda.*Ax + (1-lambda).*Bx;
Dy = lambda.*Ay + (1-lambda).*By;
D = [Dx;Dy];

end