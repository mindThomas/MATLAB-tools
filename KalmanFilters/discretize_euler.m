function [Ad, Bd] = discretize_euler(Ac, Bc, ts)
    % We discretize according to the Euler method            
    %   Ad = I + T*Ac
    %   Bd = T*Bc
    Ad = eye(size(Ac)) + ts*Ac;
    Bd = ts*Bc;    
end