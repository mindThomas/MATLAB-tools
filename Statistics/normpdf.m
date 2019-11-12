function pdf = normpdf (x, m, s)
% NORMPDF  PDF of the normal distribution
%  PDF = normpdf(X, M, S) computes the probability density
%  function (PDF) at X of the normal distribution with mean M
%  and standard deviation S.
%
%  PDF = normpdf(X) is equivalent to PDF = normpdf(X, 0, 1)

% Adapted for Matlab (R) from GNU Octave 3.0.1
% Original file: statistics/distributions/normpdf.m
% Original author: TT <Teresa.Twaroch@ci.tuwien.ac.at>

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

if (nargin ~= 1 && nargin ~= 3)
    error('normpdf: you must give one or three arguments');
end

if (nargin == 1)
    m = 0;
    s = 1;
end

if (~isscalar (m) || ~isscalar (s))
    [retval, x, m, s] = common_size (x, m, s);
    if (retval > 0)
        error ('normpdf: x, m and s must be of common size or scalars');
    end
end

sz = size (x);
pdf = zeros (sz);

if (isscalar (m) && isscalar (s))
    if (find (isinf (m) | isnan (m) | ~(s >= 0) | ~(s < Inf)))
        pdf = NaN * ones (sz);
    else
        pdf = stdnormal_pdf ((x - m) ./ s) ./ s;
    end
else
    k = find (isinf (m) | isnan (m) | ~(s >= 0) | ~(s < Inf));
    if (any (k))
        pdf(k) = NaN;
    end

    k = find (~isinf (m) & ~isnan (m) & (s >= 0) & (s < Inf));
    if (any (k))
        pdf(k) = stdnormal_pdf ((x(k) - m(k)) ./ s(k)) ./ s(k);
    end
end

pdf((s == 0) & (x == m)) = Inf;
pdf((s == 0) & ((x < m) | (x > m))) = 0;

end
