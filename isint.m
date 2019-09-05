function TF=isint(A)
%ISINT  True for integer elements.
%   Returns logical array the same size as A containing true (1)
%   where the elements of A are integers, false (0) otherwise.
%
%   TF=isint(A);

% CVS ID and authorship of this code
% CVSId = '$Id: isint.m,v 1.3 2005/02/03 16:58:34 michelich Exp $';
% CVSRevision = '$Revision: 1.3 $';
% CVSDate = '$Date: 2005/02/03 16:58:34 $';
% CVSRCSFile = '$RCSfile: isint.m,v $';

if isempty(A)
  TF=[];
elseif isnumeric(A)
  TF=imag(A)==0 & [A-floor(A)]==0;
else
  TF=logical(zeros(size(A)));
end