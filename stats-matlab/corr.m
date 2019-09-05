function retval = corr(x, y) % --*-- Unitary tests --*--
%@info:
%! @deftypefn  {Function File} {} corr (@var{x})
%! @deftypefnx {Function File} {} corr (@var{x}, @var{y})
%! Compute matrix of correlation coefficients.
%! @anchor{corr}
%!
%! If each row of @var{x} and @var{y} is an observation and each column is
%! a variable, then the @w{(@var{i}, @var{j})-th} entry of
%! @code{corr (@var{x}, @var{y})} is the correlation between the
%! @var{i}-th variable in @var{x} and the @var{j}-th variable in @var{y}.
%! @tex
%! $$
%! {\rm corr}(x,y) = {{\rm cov}(x,y) \over {\rm std}(x) {\rm std}(y)}
%! $$
%! @end tex
%! @ifnottex
%!
%! @example
%! corr (x,y) = cov (x,y) / (std (x) * std (y))
%! @end example
%!
%! @end ifnottex
%! If called with one argument, compute @code{corr (@var{x}, @var{x})},
%! the correlation between the columns of @var{x}.
%! @end deftypefn
%@eod:
%
% Notes:    - the original Octave code has been rewritten to avoid calling cov, since
%               there is a long-standing incompatiblity between Matlab's cov and Octave's cov
%               (see https://savannah.gnu.org/bugs/?40751)
%           - For compatibility with Matlab, the correlation of a constant
%               is defined as NaN, not 1
%
% Adapted for Matlab (R) from GNU Octave 4.0.1
% Original files: statistics\base\corr.m, statistics\base\cov.m, and packages\stk-2.3.4\misc\mole\corr\corr.m
% Original authors: Kurt Hornik <hornik@wu-wien.ac.at> and Julien Bect  <julien.bect@supelec.fr>

% Copyright (C) 1993-1996 Kurt Hornik
% Copyright (C) 1996-2015 John W. Eaton
% Copyright (C) 2013-2015 Julien Bect
% Copyright (C) 2016-2017 Dynare Team
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


if (nargin < 1 || nargin > 2)
    error('corr needs to be called with 1 or 2 input arguments')
end

% Input validation
if (nargin==1 && ~(isnumeric (x) || islogical (x))) || ...
           (nargin==2 && ~(isnumeric (x) || islogical (x) || isnumeric (y) || islogical (y)))
    error ('corr: X and Y must be numeric matrices or vectors');
end

if (nargin==1 && ~ismatrix(x)) || (nargin==2 && (~ismatrix(y) || ~ismatrix(x)))
    error ('corr: X and Y must be 2-D matrices or vectors');
end

if (nargin == 2)
    if size(y,1)~=size(x,1)
        error('corr: X and Y must have the same length')
    end
end

% Special case, correlation with scalar is NaN in Matlab
if isscalar(x)
    if nargin==1
        retval = NaN;
    else
        retval = NaN(size(y));
    end
    return
end

if nargin==2 && isscalar(y)
    retval = NaN(size(x'));
    return
end

n = size (x, 1);
x = x - repmat(mean(x),n,1); %demean
sx = std(x,1); %standard deviation
sx(sx==0)=NaN; %take care of constant vectors

if (nargin == 2)
    y = y - repmat (mean (y), n, 1);
    sy = std (y, 1);
    sy (sy == 0) = nan;
else
    y = x;
    sy = sx;
end

c = x'*y;
s = sx'*sy;
retval = c./(n * s);


end


%@test:1
%$ x = rand (10);
%$ cc1 = corr(x);
%$ cc2 = corr(x, x);
%$ t(1) = dassert(size(cc1),[10, 10]);
%$ t(2) = dassert(size (cc2),[10, 10]);
%$ t(3) = dassert(cc1, cc2, sqrt (eps));
%$ T = all(t);
%@eof:1

%@test:2
%$ x = [1:3]';
%$ y = [3:-1:1]';
%$ t = zeros(3,1);
%$ t(1) = dassert(corr(x, y), -1, 5*eps);
%$ t(2) = dassert(corr(x, flipud (y)), 1, 5*eps);
%$ t(3) = dassert(corr([x, y]), [1 -1; -1 1], 5*eps);
%$ T = all(t);
%@eof:2

%@test:3
%$ if ~isoctave()
%$   t = zeros(3,1);
%$   t(1) = dassert(corr(5), NaN);
%$   t(2) = dassert(corr([1 2 3],5),NaN(3,1));
%$   t(3) = dassert(corr(5,[1 2 3]),NaN(1,3));
%$ else
%$   t = 1;
%$ end
%$ T = all(t);
%@eof:3

%@test:4
%$ t = zeros(6,1);
%$ try
%$     corr()
%$     t(1) = false;
%$ catch
%$     t(1) = true;
%$ end
%$ try
%$     corr(1, 2, 3)
%$     t(2) = false;
%$ catch
%$     t(2) = true;
%$ end
%$ try
%$     corr([1; 2], ['A', 'B'])
%$     t(3) = false;
%$ catch
%$     t(3) = true;
%$ end
%$ try
%$     corr([1; 2], ['A', 'B'])
%$     t(4) = false;
%$ catch
%$     t(4) = true;
%$ end
%$ try
%$     corr(ones (2,2,2))
%$     t(5) = false;
%$ catch
%$     t(5) = true;
%$ end
%$ try
%$     corr(ones (2,2), ones (2,2,2))
%$     t(6) = false;
%$ catch
%$     t(6) = true;
%$ end
%$ T = all(t);
%@eof:4