function [Ad, Bd] = discretize_matexp(Ac, Bc, ts)
    % We discretize according to the analytic joint-matrix exponential method            
    %   M = [Ac, Bc; 0, 0]
    %   Mexp = expm(M*ts)
    %   Mexp is now = []
    M = [Ac, Bc;
         zeros(size(Bc,2), size(Ac,1)+size(Bc,2))];
    Mexp = expm(M*ts);
    % Extract discretized matrices
    Ad = Mexp(1:size(Ac,1), 1:size(Ac,2));
    Bd = Mexp(1:size(Ac,1), (size(Ac,2)+1):end);
    % This should give you the same solution as 'discretize_analytic(A,B,ts)'
end