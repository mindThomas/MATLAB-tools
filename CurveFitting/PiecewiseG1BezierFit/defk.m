function k = defk(m,n)

% Computes default knot positions for iguess.m based on a
% formula to equally disperse the knots throughout the data.

k = round(((m-1)/(n-1))*[0:n-1] + ones(1,n));