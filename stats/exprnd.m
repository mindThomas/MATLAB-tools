function rnd = exprnd(a)
%  Random samples from the exponential distribution with expectation a
%  and variance a^2.
%
%  INPUTS
%    a       [double]    m*n matrix of positive parameters
%
%  OUTPUT
%    rnd     [double]    m*n matrix, independent draws from the exponential
%                        distribution rnd(j,j) has expectation a(i,j) and
%                        variance a(i,j)^2.
%
%  ALGORITHM
%    Inverse transform sampling.
%
%  SPECIAL REQUIREMENTS
%    None.
%

% Copyright (C) 2009-2017 Dynare Team
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
if any(a(:)<0)
    disp('exprnd:: The parameter of the exponential distribution has to be positive!')
    error;
end
[m,n] = size(a);
uniform_variates = rand(m,n);
rnd = -log(uniform_variates).*a;