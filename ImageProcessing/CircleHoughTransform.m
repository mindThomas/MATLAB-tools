function [xc,yc,radius,votesNormalized] = CircleHoughTransform(edge_image)
    range_r = 4:50; % radius range
    range_theta = 0:1:359;

    votes = zeros([size(edge_image,1), size(edge_image,2), length(range_r)]);

    for (x = 1:size(edge_image,2))
        for (y = 1:size(edge_image,1))
            if (edge_image(y,x) == 1)
                for (rIdx = 1:length(range_r))
                    r = range_r(rIdx);
                    for (theta = range_theta)
                        % Construct candidate circle with center at (x,y)
                        xc = round(x + r * cos(theta * pi / 180)); % polar coordinate for center
                        yc = round(y + r * sin(theta * pi / 180));  % polar coordinate for center 
                        if (xc >= 1 && xc <= size(edge_image,2) && ...
                            yc >= 1 && yc <= size(edge_image,1))
                            votes(yc,xc,rIdx) = votes(yc,xc,rIdx) + 1; % voting
                        end
                    end
                end
            end
        end
    end

    % Find radius leading getting the most votes
    largestCount = 0;
    largestIndex = 1;
    for (idx = 1:size(votes,3))
        votes_sum = votes(:,:,idx);
        count = max(max(votes_sum));
        if (count > largestCount)
            largestCount = count;
            largestIndex = idx;
        end
    end

    radius = largestIndex + range_r(1) - 1;
    votesImage = votes(:,:,largestIndex);
    votesNormalized = votesImage / max(max(votesImage));

    [maxval,idx] = max(votesNormalized(:));
    [yc,xc] = ind2sub(size(votesNormalized), idx);
end