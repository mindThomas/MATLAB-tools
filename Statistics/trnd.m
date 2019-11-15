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
%% @deftypefn  {} {} trnd (@var{n})
%% @deftypefnx {} {} trnd (@var{n}, @var{r})
%% @deftypefnx {} {} trnd (@var{n}, @var{r}, @var{c}, @dots{})
%% @deftypefnx {} {} trnd (@var{n}, [@var{sz}])
%% Return a matrix of random samples from the t (Student) distribution with
%% @var{n} degrees of freedom.
%%
%% When called with a single size argument, return a square matrix with
%% the dimension specified.  When called with more than one scalar argument the
%% first two arguments are taken as the number of rows and columns and any
%% further arguments specify additional matrix dimensions.  The size may also
%% be specified with a vector of dimensions @var{sz}.
%%
%% If no size arguments are given then the result matrix is the size of
%% @var{n}.
%% @end deftypefn

%% Author: KH <Kurt.Hornik@wu-wien.ac.at>
%% Description: Random deviates from the t distribution

function rnd = trnd (n, varargin)

  if (nargin < 1)
    error("trnd: Incorrect usage: rnd = trnd(n, varargin)");
  end

  if (nargin == 1)
    sz = size (n);
  elseif (nargin == 2)
    if (isscalar (varargin{1}) && varargin{1} >= 0)
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0))
      sz = varargin{1};
    else
      error ("trnd: dimension vector must be row vector of non-negative integers");
    end
  elseif (nargin > 2)
    if (any (cellfun (@(x) (~ isscalar (x) || x < 0), varargin)))
      error ("trnd: dimensions must be non-negative integers");
    end
    sz = [varargin{:}];
  end

  if (~ isscalar (n) && ~ isequal (size (n), sz))
    error ("trnd: N must be scalar or of size SZ");
  end

  if (~isreal (n))
    error ("trnd: N must not be complex");
  end

  if (isa (n, "single"))
    cls = "single";
  else
    cls = "double";
  end

  if (isscalar (n))
    if ((n > 0) && (n < Inf))
      gammaRandom = gamrnd(ones(sz)*n/2, ones(sz)*1);
      gammaRandom = reshape(gammaRandom, sz);
        
      rnd = randn (sz, cls) ./ sqrt (2*gammaRandom ./ n);
    else
      rnd = NaN (sz, cls);
    end
  else
    rnd = NaN (sz, cls);

    k = (n > 0) & (n < Inf);
    rnd(k) = randn (sum (k(:)), 1, cls) ...
             ./ sqrt (2*randg (n(k)/2, cls) ./ n(k));
  end

end