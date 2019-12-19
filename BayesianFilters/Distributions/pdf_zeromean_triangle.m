function p = pdf_zeromean_triangle(x, var)
    % Evaluate the probability (PDF) of a zero-centered (standard) Triangle
    % distribution with variance, var
    % See Table 5.2 in Probabilistic Robotics    
    p = max(0,  1/sqrt(6*var) - abs(a)/(6*var)   );
end