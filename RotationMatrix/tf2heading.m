function heading = tf2heading(tf, varargin)    
    rotm = tf2rotm(tf);
    heading = rotm2heading(rotm);
end
