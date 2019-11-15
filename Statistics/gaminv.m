function inv = gaminv (x, a, b)
% GAMINV  Quantile function of the Gamma distribution
%  INV = gaminv(X, A, B) computes, for each component of X, the
%  quantile (the inverse of the CDF) at X of the Gamma distribution
%  with parameters A and B (i.e. mean of the distribution is A*B
%  and variance is A*B^2).

% Adapted for Matlab (R) from GNU Octave 3.0.1
% Original file: statistics/distributions/gaminv.m
% Original author: KH <Kurt.Hornik@wu-wien.ac.at>

% Copyright (C) 1995, 1996, 1997, 2005, 2006, 2007 Kurt Hornik
% Copyright (C) 2008-2009 Dynare Team
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
    error ('gaminv: you must give three arguments');
end

if (~isscalar (a) || ~isscalar(b))
    [retval, x, a, b] = common_size (x, a, b);
    if (retval > 0)
        error ('gaminv: x, a and b must be of common size or scalars');
    end
end

sz = size (x);
inv = zeros (sz);

k = find ((x < 0) | (x > 1) | isnan (x) | ~(a > 0) | ~(b > 0));
if (any (k))
    inv (k) = NaN;
end

k = find ((x == 1) & (a > 0) & (b > 0));
if (any (k))
    inv (k) = Inf;
end

k = find ((x > 0) & (x < 1) & (a > 0) & (b > 0));
if (any (k))
    if (~isscalar(a) || ~isscalar(b))
        a = a (k);
        b = b (k);
        y = a .* b;
    else
        y = a * b * ones (size (k));
    end
    x = x (k);
    l = find (x < eps);
    if (any (l))
        y(l) = sqrt (eps) * ones (length (l), 1);
    end

    y_old = y;
    for i = 1 : 100
        h     = (gamcdf (y_old, a, b) - x) ./ gampdf (y_old, a, b);
        y_new = y_old - h;
        ind   = find (y_new <= eps);
        if (any (ind))
            y_new (ind) = y_old (ind) / 10;
            h = y_old - y_new;
        end
        if (max (abs (h)) < sqrt (eps))
            break
        end
        y_old = y_new;
    end

    inv (k) = y_new;
end

end
