function p = wblcdf(x, scale, shape)   % --*-- Unitary tests --*--

% Cumulative distribution function for the Weibull distribution.
%
% INPUTS
% - x     [double] Positive real scalar.
% - scale [double] Positive hyperparameter.
% - shape [double] Positive hyperparameter.
%
% OUTPUTS
% - p     [double] Positive scalar between

% Copyright (C) 2015-2017 Dynare Team
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


% Check input arguments.

if nargin<3
    error('Three input arguments required!')
end

if ~isnumeric(x) || ~isscalar(x) || ~isreal(x)
    error('First input argument must be a real scalar!')
end

if ~isnumeric(scale) || ~isscalar(scale) || ~isreal(scale) || scale<=0
    error('Second input argument must be a real positive scalar (scale parameter of the Weibull distribution)!')
end

if ~isnumeric(shape) || ~isscalar(shape) || ~isreal(shape) || shape<=0
    error('Third input argument must be a real positive scalar (shape parameter of the Weibull distribution)!')
end

% Filter trivial polar cases.

if x<=0
    p = 0;
    return
end

if isinf(x)
    p = 1;
    return
end

% Evaluate the CDF.

p = 1-exp(-(x/scale)^shape);

%@test:1
%$ try
%$     p = wblcdf(-1, .5, .1);
%$     t(1) = true;
%$ catch
%$     t(1) = false;
%$ end
%$
%$ % Check the results
%$ if t(1)
%$     t(2) = isequal(p, 0);
%$ end
%$ T = all(t);
%@eof:1

%@test:2
%$ try
%$     p = wblcdf(Inf, .5, .1);
%$     t(1) = true;
%$ catch
%$     t(1) = false;
%$ end
%$
%$ % Check the results
%$ if t(1)
%$     t(2) = isequal(p, 1);
%$ end
%$ T = all(t);
%@eof:2

%@test:3
%$ % Set the hyperparameters of a Weibull definition.
%$ scale = .5;
%$ shape = 1.5;
%$
%$ % Compute the median of the weibull distribution.
%$ m = scale*log(2)^(1/shape);
%$
%$ try
%$     p = wblcdf(m, scale, shape);
%$     t(1) = true;
%$ catch
%$     t(1) = false;
%$ end
%$
%$ % Check the results
%$ if t(1)
%$     t(2) = abs(p-.5)<1e-12;
%$ end
%$ T = all(t);
%@eof:3

%@test:4
%$ % Consistency check between wblinv and wblcdf.
%$
%$ % Set the hyperparameters of a Weibull definition.
%$ scale = .5;
%$ shape = 1.5;
%$
%$ % Compute quatiles of the weibull distribution.
%$ q = 0:.05:1;
%$ m = zeros(size(q));
%$ p = zeros(size(q));
%$ for i=1:length(q)
%$     m(i) = wblinv(q(i), scale, shape);
%$ end
%$
%$ try
%$     for i=1:length(q)
%$         p(i) = wblcdf(m(i), scale, shape);
%$     end
%$     t(1) = true;
%$ catch
%$     t(1) = false;
%$ end
%$
%$ % Check the results
%$ if t(1)
%$     for i=1:length(q)
%$         t(i+1) = abs(p(i)-q(i))<1e-12;
%$     end
%$ end
%$ T = all(t);
%@eof:4
