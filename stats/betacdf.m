function cdf = betacdf (x, a, b)
% BETACDF  CDF of the Beta distribution
%  CDF = betacdf(X, A, B) computes, for each element of X, the CDF
%  at X of the beta distribution with parameters A and B (i.e.
%  mean of the distribution is A/(A+B) and variance is
%  A*B/(A+B)^2/(A+B+1) ).

% Adapted for Matlab (R) from GNU Octave 3.0.1
% Original file: statistics/distributions/betacdf.m
% Original author: KH <Kurt.Hornik@wu-wien.ac.at>

% Copyright (C) 1995, 1996, 1997, 2005, 2006, 2007 Kurt Hornik
% Copyright (C) 2008-2017 Dynare Team
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
    error ('betacdf: you should provide three arguments');
end

if (~isscalar (a) || ~isscalar(b))
    [retval, x, a, b] = common_size (x, a, b);
    if (retval > 0)
        error ('betacdf: x, a and b must be of common size or scalar');
    end
end

sz = size(x);
cdf = zeros (sz);

k = find (~(a > 0) | ~(b > 0) | isnan (x));
if (any (k))
    cdf (k) = NaN;
end

k = find ((x >= 1) & (a > 0) & (b > 0));
if (any (k))
    cdf (k) = 1;
end

k = find ((x > 0) & (x < 1) & (a > 0) & (b > 0));
if (any (k))
    if (isscalar (a) && isscalar(b))
        cdf (k) = betainc (x(k), a, b);
    else
        cdf (k) = betainc (x(k), a(k), b(k));
    end
end

end
