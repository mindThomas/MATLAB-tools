function eigs = ordeig(t)
% function eval = ordeig(t)
% Computes the eigenvalues of a quasi-triangular matrix
%
% INPUTS
%    t:              quasi-triangular matrix
%
% OUTPUTS
%    eigs:           eigenvalues
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2003-2017 Dynare Team
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

n = size(t,2);
eigs = zeros(n,1);
i = 1;
while i <= n
    if i == n
        eigs(n) = t(n,n);
        break
    elseif t(i+1,i) == 0
        eigs(i) = t(i,i);
        i = i+1;
    else
        k = i:i+1;
        eigs(k) = eig(t(k,k));
        i = i+2;
    end
end
