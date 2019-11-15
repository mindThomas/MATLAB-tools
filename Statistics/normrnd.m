% Copyright (C) 2012 Rik Wehbring
% Copyright (C) 1995-2016 Kurt Hornik
%
% This program is free software: you can redistribute it and/or
% modify it under the terms of the GNU General Public License as
% published by the Free Software Foundation, either version 3 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; see the file COPYING.  If not, see
% <http://www.gnu.org/licenses/>.

% -*- texinfo -*-
% @deftypefn  {} {} normrnd (@var{mu}, @var{sigma})
% @deftypefnx {} {} normrnd (@var{mu}, @var{sigma}, @var{r})
% @deftypefnx {} {} normrnd (@var{mu}, @var{sigma}, @var{r}, @var{c}, @dots{})
% @deftypefnx {} {} normrnd (@var{mu}, @var{sigma}, [@var{sz}])
% Return a matrix of random samples from the normal distribution with
% parameters mean @var{mu} and standard deviation @var{sigma}.
%
% When called with a single size argument, return a square matrix with
% the dimension specified.  When called with more than one scalar argument the
% first two arguments are taken as the number of rows and columns and any
% further arguments specify additional matrix dimensions.  The size may also
% be specified with a vector of dimensions @var{sz}.
%
% If no size arguments are given then the result matrix is the common size of
% @var{mu} and @var{sigma}.
% @end deftypefn

% Author: KH <Kurt.Hornik@wu-wien.ac.at>
% Description: Random deviates from the normal distribution

function rnd = normrnd (mu, sigma, varargin)
    if (nargin < 2)
        error('Incorrect usage');  
    end

    if (~isscalar (mu) || ~isscalar (sigma))
        [retval, mu, sigma] = common_size (mu, sigma);
        if (retval > 0)
            error ("normrnd: MU and SIGMA must be of common size or scalars");
        end
    end

    if (iscomplex (mu) || iscomplex (sigma))
        error ("normrnd: MU and SIGMA must not be complex");
    end

    if (nargin == 2)
        sz = size (mu);
    elseif (nargin == 3)
        if (isscalar (varargin{1}) && varargin{1} >= 0)
            sz = [varargin{1}, varargin{1}];
        elseif (isrow (varargin{1}) && all (varargin{1} >= 0))
            sz = varargin{1};
        else
            error ("normrnd: dimension vector must be row vector of non-negative integers");
        end
    elseif (nargin > 3)
        if (any (cellfun (@(x) (~isscalar (x) || x < 0), varargin)))
            error ("normrnd: dimensions must be non-negative integers");
        end
        sz = [varargin{:}];
    end

    if (~isscalar (mu) && ~isequal (size (mu), sz))
        error ("normrnd: MU and SIGMA must be scalar or of size SZ");
    end

    if (isa (mu, "single") || isa (sigma, "single"))
        cls = "single";
    else
        cls = "double";
    end

    if (isscalar (mu) && isscalar (sigma))
        if (isfinite (mu) && (sigma >= 0) && (sigma < Inf))
            rnd = mu + sigma * randn (sz, cls);
        else
            rnd = NaN (sz, cls);
        end
    else
        rnd = mu + sigma .* randn (sz, cls);
        k = ~isfinite (mu) | ~(sigma >= 0) | ~(sigma < Inf);
        rnd(k) = NaN;
    end

end