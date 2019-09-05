function rnd = betarnd(a, b)
% This function produces independent random variates from the Beta distribution.
%
%  INPUTS
%    a       [double]    n*1 vector of positive parameters.
%    b       [double]    n*1 vector of positive parameters.
%
%  OUTPUT
%    rnd     [double]    n*1 vector of independent variates from the beta(a,b) distribution.
%                        rnd(i) is beta distributed with variance a(i)/(a(i)+b(i)) and
%                        variance a(i)b(i)/(a(i)+b(i))^2/(a(i)+b(i)+1).
%
%  ALGORITHMS
%    Described and Devroye (1986, chapter 9).

% Copyright (C) 2008-2017 Dynare Team
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
    error('betarnd: you must give two arguments');
end

if (any(a<0)) || (any(b<0)) || (any(a==Inf)) || (any(b==Inf))
    error('betarnd:: Input arguments must be finite and positive!');
end

[ma,na] = size(a);
[mb,nb] = size(b);

if ma~=mb || na~=nb
    error('betarnd:: Input arguments must have the same size!');
end

if na~=1
    error('betarnd:: Input arguments must be column vectors');
end

x = gamrnd(a,ones(ma,1));
rnd = x./(x+gamrnd(b, ones(mb,1)));