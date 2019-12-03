function x = rand_zeromean_triangle(var)
    % See Table 5.4 from Probabilistic Robotics
    stddev = sqrt(var);
    r = 2*stddev*rand(2,1) - stddev;
    x = sqrt(6)/2 * sum(r);
end