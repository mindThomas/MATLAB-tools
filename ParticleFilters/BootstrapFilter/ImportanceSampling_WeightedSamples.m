function [x, weights] = ImportanceSampling_WeightedSamples(target_pdf, proposal_distribution, n_realizations)
    % See "7.2.1 Monte Carlo (MC) and Importance Sampling (IS)" from Course ChM015x
    % 
    % Use Importance Sampling to draw approximate realizations from a 
    % target distribution from which it is otherwise hard to draw true samples
    % target_pdf is a PDF function of type @(x) ...
    % proposal_distribution is a distribution object which includes the
    % function .draw()
    %
    % Note that it is important that the support of the proposal
    % distribution includes the support of the target distribution
    % 
    % Note that importance sampling does not itself draw the final
    % approximating samples, but instead weights each of the drawn samples
    % according to the distribution ratio.
    % Afterwards Resampling can be applied to get the actual approximating
    % samples   
    
    x = [];
    weights = [];
    
    % Draw realizations from the proposal distribution and compute weight
    % from the distribution ratio
    for (i = 1:n_realizations)
        sample = proposal_distribution.draw();
        x(:,end+1) = sample;
        weights(end+1) = target_pdf(sample) / proposal_distribution.pdf(sample);
    end
    
    % Normalize weights
    weights = weights / sum(weights);
end