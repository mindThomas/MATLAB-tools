function [Ad, Bd] = discretize_analytic(Ac, Bc, ts)
    % We discretize according to the analytic method with zero order hold input
    %   Ad = exp(Ac*ts)     [ = expm(Ac*ts) ]
    %   Bd = int_0^ts exp(A*tau) dtau * Bc
    % The integral can be computed from the Taylor expansion of exp(A*t)
    %   exp(A*t) = I + A*t + A^2*t^2/2 + A^3*t^3/3! + A^4*t^4/4! + ...
    % such that
    %   int_0^ts exp(A*tau) dtau = I*ts + A*ts^2/2 + A^2*ts^3/3! + ...
    Ad = eye(size(Ac));
    Bd = eye(size(Ac)) * ts;
    
    A_mul = Ac;
    ts_mul = ts;
    i_factorial = 1;    
    
    for (i = 2:50)        
        Ad = Ad + A_mul * ts_mul / i_factorial;
        i_factorial = i_factorial * i;
        ts_mul = ts_mul * ts;
        Bd = Bd + A_mul * ts_mul / i_factorial;
        A_mul = A_mul * Ac;
    end
    if (~isempty(Bc))
        Bd = Bd * Bc;
    else
        Bd = [];
    end
end