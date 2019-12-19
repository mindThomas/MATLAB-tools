function x = rand_zeromean_gaussian_efficient(var)
    % See Table 5.4 from Probabilistic Robotics
    stddev = sqrt(var);
    r = 2*stddev*rand(12,1) - stddev;
    x = 1/2 * sum(r);
end