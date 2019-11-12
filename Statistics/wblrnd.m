function rnd = wblrnd(a, b)
% This function produces independent random variates from the Weibull distribution.
%
%  INPUTS
%    a       [double]    m*n matrix of positive parameters (scale).
%    b       [double]    m*n matrix of positive parameters (shape).
%
%  OUTPUT
%    rnd     [double]    m*n matrix of independent variates from the beta(a,b) distribution.

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

if (nargin ~= 2)
    error('Two input arguments required!');
end

if (any(a(:)<0)) || (any(b(:)<0)) || (any(a(:)==Inf)) || (any(b(:)==Inf))
    error('Input arguments must be finite and positive!');
end

[ma,na] = size(a);
[mb,nb] = size(b);

if ma~=mb || na~=nb
    error('Input arguments must have the same size!');
end

rnd = a.*(-log(rand(ma, na))).^(1./b);