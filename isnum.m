function TF=isnum(A)
%ISNUM  True for numeric elements.
%   Returns logical array the same size as A containing true (1)
%   where the elements of A are integers, false (0) otherwise.
%
%   TF=isnum(A);

TF=logical(zeros(size(A)));
for (i = 1:length(A))
    TF(i) = (isnumeric(A(i)) && (imag(A(i))) == 0);
end