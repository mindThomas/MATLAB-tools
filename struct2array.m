function a = struct2array(s)

% INPUTS
% - s  [struct]  with N fields, field i contains a n_i*m_i array of doubles.
%
% OUPUTS
% - a  [double]  column vector with sum(n_i*m_i, i=1,...,N) elements.

% Copyright (C) 2017 Dynare Team
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

assert(isstruct(s), 'struct2array:: Argument has to be a structure!')

c = cellfun(@vec, struct2cell(s), 'UniformOutput', false);
a = vertcat(c{:});