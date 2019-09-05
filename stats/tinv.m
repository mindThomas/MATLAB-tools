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
%% @deftypefn {} {} tinv (@var{x}, @var{n})
%% For each element of @var{x}, compute the quantile (the inverse of the CDF)
%% at @var{x} of the t (Student) distribution with @var{n}
%% degrees of freedom.
%%
%% This function is analogous to looking in a table for the t-value of a
%% single-tailed distribution.
%% @end deftypefn

%% For very large n, the "correct" formula does not really work well,
%% and the quantiles of the standard normal distribution are used
%% directly.

%% Author: KH <Kurt.Hornik@wu-wien.ac.at>
%% Description: Quantile function of the t distribution

function inv = tinv (x, n)

  if (nargin ~= 2)
    error("tinv: Incorrect usage: inv = tinv(x, n)");
  end

  if (~isreal (x) || ~isreal (n))
    error ("tinv: X and N must not be complex");
  end

  if (isa (x, "single") || isa (n, "single"))
    inv = NaN (size (x), "single");
  else
    inv = NaN (size (x));
  end

  k = (x == 0) & (n > 0);
  inv(k) = -Inf;

  k = (x == 1) & (n > 0);
  inv(k) = Inf;

  if (isscalar (n))
    k = (x > 0) & (x < 1);
    if ((n > 0) && (n < 10000))
      inv(k) = (sign (x(k) - 1/2) ...
                .* sqrt (n * (1 ./ betainv (2*min (x(k), 1 - x(k)), ...
                                            n/2, 1/2) - 1)));
    elseif (n >= 10000)
      %% For large n, use the quantiles of the standard normal
      inv(k) = stdnormal_inv (x(k));
    end
  else
    k = (x > 0) & (x < 1) & (n > 0) & (n < 10000);
    inv(k) = (sign (x(k) - 1/2) ...
              .* sqrt (n(k) .* (1 ./ betainv (2*min (x(k), 1 - x(k)), ...
                                              n(k)/2, 1/2) - 1)));

    %% For large n, use the quantiles of the standard normal
    k = (x > 0) & (x < 1) & (n >= 10000);
    inv(k) = stdnormal_inv (x(k));
  end

end
