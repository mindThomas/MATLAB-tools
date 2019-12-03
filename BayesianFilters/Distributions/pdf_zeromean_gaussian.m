function p = pdf_zeromean_gaussian(x, var)
    % Evaluate the probability (PDF) of a zero-centered (standard) Gaussian
    % distribution with variance, var
    % See Table 5.2 in Probabilistic Robotics
    p = 1/(2*pi*var) * exp(-1/(2*var) * x^2);
end