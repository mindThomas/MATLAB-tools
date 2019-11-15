function pdf = stdnormal_pdf (x)
% STDNORMAL_PDF  PDF of the standard normal distribution
%  PDF = stdnormal_pdf(X)
%  For each element of X, compute the PDF of the standard normal
%  distribution at X.

% Adapted for Matlab (R) from GNU Octave 3.0.1
% Original file: statistics/distributions/stdnormal_pdf.m
% Original author: TT <Teresa.Twaroch@ci.tuwien.ac.at>

% Copyright (C) 1995, 1996, 1997, 1998, 2000, 2002, 2004, 2005, 2006,
%               2007 Kurt Hornik
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

if (nargin ~= 1)
    error('stdnormal_pdf: you should provide one argument');
end

sz = size(x);
pdf = zeros (sz);

k = find (isnan (x));
if (any (k))
    pdf(k) = NaN;
end

k = find (~isinf (x));
if (any (k))
    pdf (k) = (2 * pi)^(- 1/2) * exp (- x(k) .^ 2 / 2);
end

end
