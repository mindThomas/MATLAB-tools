function tok = strsplit(str, delimiters)

% Splits a string into multiple terms.
%
% INPUTS
% - str        [char]                String to be splitted.
% - delimiters [char, cell(char)]    Delimiters.
%
% OUTPUTS
% - tok        [cell(char)]          Terms.

% Copyright © 2018 DynareTeam
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

remove_empty = true;
remove_numbers = false;

% Check first input arguments
assert(ischar(str) && ndims(str)==2 && size(str,1)<=1, 'The first arugment has to be a row char array!');

% Set default value for second input arguments
if nargin<2
    delimiters = {' '};
end

% If second input argument is a char transform it into a sigleton cell of char
if nargin>1
    if ischar(delimiters)
        assert(ndims(delimiters)==2 && size(delimiters,1)==1, 'The second input argument has to be be a char string!');
        delimiters = {delimiters};
    end
end

% Check that `delimiters` is a one dimensional cell
assert(ndim(delimiters)<=1, 'The second input argument has to be a one dimensional cell array!')

% Check that `delimiters` is a cell of row char arrays
assert(all(cellfun(@ischar, delimiters)) && all(cellfun(@rows, delimiters)==1), 'The second input argument has to be a cell of row char arrays!')

% If space is one of the delimiters obtain the index in `delimiters`
idspace = strmatch(' ', delimiters);

% Get the number of delimiters
n = length(delimiters);

% Remove unnecessary spaces
delimiters(setdiff(1:n, idspace)) = strtrim(delimiters(setdiff(1:n, idspace)));

% Join all the delimiters (strjoin is not available with matlab version less than R2013a)
if n>1
    delimiter = '';
    for i=1:n
        if isspace(delimiters{i})
            delimiter = horzcat(delimiter, '\s');
        else
            delimiter = horzcat(delimiter, delimiters{i});
        end
        delimiter = horzcat(delimiter,'|');
    end
    delimiter = horzcat(delimiter, '\W');
else
    delimiter = delimiters{1};
end

% Get tokens
tok = regexp(str, delimiter, 'split');

if remove_empty
    % Remove empty tokens
    tok = tok(find(~cellfun(@isempty, tok)));
end

if remove_numbers
    % Remove numbers
    tok = tok(find(~cellfun(@ischarnum, tok)));
end