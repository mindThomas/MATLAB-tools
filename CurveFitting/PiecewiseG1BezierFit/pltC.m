function graf = pltC(C,Q,P)

% function graf = pltC(C,Q,P). This function takes as input:
% control points,C, data points,Q; and knot points, P. The
% control points are used to calculate the points of the
% approximating cubic Bezier curves. The control, data, and
% knot points are then plotted along with the curve.
% This was written by M. R. Holmes.

[s,t] = size(C);
x = [0:0.025:1]; % Defines the interval for the polynomial.
[a,b] = size(x);

W = [ ]; % Loop to construct the Bezier curve.
for j = 1:3:t-3
	Y = zeros(2,b);
	M = [berny(3,0,x)' berny(3,1,x)' berny(3,2,x)' berny(3,3,x)'];
	Y = Y + C(:,j:j+3) * M';
	W = [W Y];
end

plot( W(1,:) , W(2,:)); hold on;
plot( C(1,:) , C(2,:));
plot( Q(1,:) , Q(2,:), '+');
plot( P(1,:) , P(2,:), 'x');
plot( C(1,:) , C(2,:), 'o');

graf = gcf;