function inv = chi2inv (x, n)
% CHI2INV  Quantile function of the chi-square distribution
%  INV = chi2inv(X, N) computes, for each element of X, the
%  quantile (the inverse of the CDF) at X of the chi-square
%  distribution with N degrees of freedom.

% Adapted for Matlab (R) from GNU Octave 3.0.1
% Original file: statistics/distributions/chi2inv.m
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

if (nargin ~= 2)
    error ('chi2inv: you must give two arguments');
end

if (~isscalar (n))
    [retval, x, n] = common_size (x, n);
    if (retval > 0)
        error ('chi2inv: x and n must be of common size or scalar');
    end
end

inv = gaminv (x, n / 2, 2);

end
