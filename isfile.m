function a = isfile(b)

%@info:
%! @deftypefn {Function File} {@var{a} =} isfile (@var{b})
%! @anchor{isfile}
%! @sp 1
%! Test if @var{b} is a file.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @var
%! @item b
%! A matlab/octave string.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @var
%! @item a
%! Integer scalar, equal to 1 if @var{b} is a file, zero otherwise.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%! @sp 2
%! @end deftypefn
%@eod:

% Copyright (C) 2012-2017 Dynare Team
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

stringarrayflag = false;
cellofstringflag = false;
n = 1;
a = false;

if ~isoctave() && ~matlab_ver_less_than('9.1') && isstring(b) && length(b)>1 && isvector(b)
    n = length(b);
    stringarrayflag = true;
    a = false(size(b));
end

if iscell(b) && length(b)>1 && isvector(b)
    if all(cellfun(@ischar, b))
        n = length(b);
        cellofstringflag = true;
        a = false(size(b));
    else
        error('Wrong input argument type!')
    end
end

for i=1:n
    if stringarrayflag
        d = b(i);
    elseif cellofstringflag
        d = b{i};
    elseif ischar(b) && size(b, 1)==1 
        d = b;
    else
        error('Wrong input argument type!')
    end
    [base, ext] = strtok(d, '.');
    if isempty(ext)
        % File has no extension.
        [status, c] = fileattrib(d);
        if status
            a(i) = ~c.directory;
        end
    else
        a(i) = isequal(exist(d, 'file'), 2);
    end
end