function [Ad, Bd] = discretize_bilinear(Ac, Bc, ts)    
    % We discretize by applying the bilinear transform   
    % See https://dsp.stackexchange.com/questions/45042/bilinear-transformation-of-continuous-time-state-space-system
	IT = (2/ts) * eye(size(Ac));
    invITA = inv(IT-Ac);
    Ad = invITA * (IT+Ac);
    iab = invITA * Bc;
    tk = 2 / sqrt(ts);
    Bd = 2*iab; % not sure if this is correct
    %Cd = tk*(Cc/(IT-Ac));
    %Dd = Dc + (Cc*iab);
end