function pdf = betapdf (x, a, b)
% BETAPDF  PDF of the Beta distribution
%  PDF = betapdf(X, A, B) computes, for each element of X, the PDF
%  at X of the beta distribution with parameters A and B (i.e.
%  mean of the distribution is A/(A+B) and variance is
%  A*B/(A+B)^2/(A+B+1) ).

% Adapted for Matlab (R) from GNU Octave 3.0.1
% Original file: statistics/distributions/betapdf.m
% Original author: KH <Kurt.Hornik@wu-wien.ac.at>
% Modified by Michel Juillard <michel.juillard@mjui.fr> for large values of a and b

% Copyright (C) 1995, 1996, 1997, 2005, 2006, 2007 Kurt Hornik
% Copyright (C) 2008-2011 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

if (nargin ~= 3)
    error ('betapdf: you must give three arguments');
end

if (~isscalar (a) || ~isscalar(b))
    [retval, x, a, b] = common_size (x, a, b);
    if (retval > 0)
        error ('betapdf: x, a and b must be of common size or scalar');
    end
end

sz = size (x);
pdf = zeros (sz);

k = find (~(a > 0) | ~(b > 0) | isnan (x));
if (any (k))
    pdf (k) = NaN;
end

k = find ((x > 0) & (x < 1) & (a > 0) & (b > 0));
if (any (k))
    if (isscalar(a) && isscalar(b))
        pdf(k) = exp ((a - 1) .* log (x(k)) ...
                      + (b - 1) .* log (1 - x(k))-gammaln(a)-gammaln(b)+gammaln(a+b));
    else
        pdf(k) = exp ((a(k) - 1) .* log (x(k)) ...
                      + (b(k) - 1) .* log (1 - x(k))-gammaln(a(k))-gammaln(b(k))+gammaln(a(k)+b(k)));
    end
end

end
