% Take scalar power n of a quaternion
% https://math.stackexchange.com/a/939288
function q_out = quatpow(q, n)
    q_out = quatexp( n * quatlogn(q) );
end