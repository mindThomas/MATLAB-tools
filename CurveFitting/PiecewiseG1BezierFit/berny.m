function val = berny(n, i, x)

% This function is a non-recursive formula for cubic Bernstein
% Polynomials which form the basis for the cubic Bezier curves
% which are used in the supporting programs. The inputs are the
% degree of the polynomial, n; the particular curve that is ass-
% igned a value of zero up to and including the degree, ii and
% the points between [0,1) to be evaluated, x. The output is
% points on the curve. It was written by M. R. Holmes.

ni = [1 3 3 1]; 
m = size(x);
if n < i
	val = zeros (m) ;
elseif i < 0
	val = zeros (m) ;
elseif ((n == 0) && (i == 0))
	val = 1;
else
	val = ni(i+1) * (x.^i) .* ((ones(m) - x) .^(n-i));
end