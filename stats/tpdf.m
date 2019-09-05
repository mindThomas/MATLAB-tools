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
%% @deftypefn {} {} tpdf (@var{x}, @var{n})
%% For each element of @var{x}, compute the probability density function (PDF)
%% at @var{x} of the @var{t} (Student) distribution with
%% @var{n} degrees of freedom.
%% @end deftypefn

%% Author: KH <Kurt.Hornik@wu-wien.ac.at>
%% Description: PDF of the t distribution

function pdf = tpdf (x, n)

  if (nargin ~= 2)
    error("tpdf: Incorrect usage: pdf = tpdf(x, n)");
  end

  if (~ isscalar (n))
    [retval, x, n] = common_size (x, n);
    if (retval > 0)
      error ("tpdf: X and N must be of common size or scalars");
    end
  end

  if (~isreal (x) || ~isreal (n))
    error ("tpdf: X and N must not be complex");
  end

  if (isa (x, "single") || isa (n, "single"))
    pdf = zeros (size (x), "single");
  else
    pdf = zeros (size (x));
  end

  k = isnan (x) | ~(n > 0) | ~(n < Inf);
  pdf(k) = NaN;

  k = isfinite (x) & (n > 0) & (n < Inf);
  if (isscalar (n))
    pdf(k) = (exp (- (n + 1) * log (1 + x(k) .^ 2 / n)/2) ...
              / (sqrt (n) * beta (n/2, 1/2)));
  else
    pdf(k) = (exp (- (n(k) + 1) .* log (1 + x(k) .^ 2 ./ n(k))/2) ...
              ./ (sqrt (n(k)) .* beta (n(k)/2, 1/2)));
  end

end