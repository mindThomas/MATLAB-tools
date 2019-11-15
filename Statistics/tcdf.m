%% Copyright (C) 2013-2017 Julien Bect
%% Copyright (C) 2012 Rik Wehbring
%% Copyright (C) 1995-2016 Kurt Hornik
%%
%% This program is free software: you can redistribute it and/or
%% modify it under the terms of the GNU General Public License as
%% published by the Free Software Foundation, either version 3 of the
%% License, or (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful, but
%% WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%% General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; see the file COPYING.  If not, see
%% <http://www.gnu.org/licenses/>.

%% -*- texinfo -*-
%% @deftypefn {} {} tcdf (@var{x}, @var{n})
%% For each element of @var{x}, compute the cumulative distribution function
%% (CDF) at @var{x} of the t (Student) distribution with
%% @var{n} degrees of freedom.
%% @end deftypefn

%% Author: KH <Kurt.Hornik@wu-wien.ac.at>
%% Description: CDF of the t distribution

function cdf = tcdf (x, n)

  if (nargin ~= 2)
    error("tcdf: Incorrect usage: cdf = tcdf(x, n)");
  end

  if (~ isscalar (n))
    [retval, x, n] = common_size (x, n);
    if (retval > 0)
      error ("tcdf: X and N must be of common size or scalars");
    end
  end

  if (~isreal (x) || ~isreal (n))
    error ("tcdf: X and N must not be complex");
  end

  if (isa (x, "single") || isa (n, "single"))
    cdf = zeros (size (x), "single");
  else
    cdf = zeros (size (x));
  end

  k = ~ isinf (x) & (n > 0);

  xx = x .^ 2;
  x_big_abs = (xx > n);

  %% deal with the case "abs(x) big"
  kk = k & x_big_abs;
  if (isscalar (n))
    cdf(kk) = betainc (n ./ (n + xx(kk)), n/2, 1/2) / 2;
  else
    cdf(kk) = betainc (n(kk) ./ (n(kk) + xx(kk)), n(kk)/2, 1/2) / 2;
  end

  %% deal with the case "abs(x) small"
  kk = k & ~ x_big_abs;
  if (isscalar (n))
    cdf(kk) = 0.5 * (1 - betainc (xx(kk) ./ (n + xx(kk)), 1/2, n/2));
  else
    cdf(kk) = 0.5 * (1 - betainc (xx(kk) ./ (n(kk) + xx(kk)), 1/2, n(kk)/2));
  end

  k = k .* (x > 0);
  if (any (k(:)))
    cdf(k) = 1 - cdf(k);
  end

  k = isnan (x) | ~(n > 0);
  cdf(k) = NaN;

  k = (x == Inf) & (n > 0);
  cdf(k) = 1;

end