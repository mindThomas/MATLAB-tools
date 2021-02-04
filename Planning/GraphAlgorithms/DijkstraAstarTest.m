raw_track = imread('track.png');
track = raw_track(:,:,1) > 0;

figure(1);
imshow(255 * track);

%% Define Dijkstra neighbour cost
neighbor_cost_matrix = inf * ones(5,5);
c = (((size(neighbor_cost_matrix)-[1,1])/2)+[1,1]);
for i = 1:size(neighbor_cost_matrix,1)
    for j = 1:size(neighbor_cost_matrix,2)
        if ~((i == c(1) && j == c(2)) || ...
             (abs(i-c(1)) > 1 && abs(i-c(1)) == abs(j-c(2))) || ...
             (abs(i-c(1)) == 0 && abs(j-c(2)) > 1) || ...
             (abs(j-c(2)) == 0 && abs(i-c(1)) > 1) ...
             )
            neighbor_cost_matrix(i, j) = sqrt( (i-c(1))^2 + (j-c(2))^2 );
        end
    end
end

% 5x5 neighbor grid, using Euclidean distance as cost
% neighbor_cost_matrix = [inf, sqrt(1^2 + 2^2), sqrt(0^2 + 2^2), sqrt(1^2 + 2^2), inf
%                         sqrt(2^2 + 1^2), sqrt(1^2 + 1^2), 1, sqrt(1^2 + 1^2), sqrt(2^2 + 1^2)
%                         inf, 1, inf, 1, inf
%                         sqrt(2^2 + 1^2), sqrt(1^2 + 1^2), 1, sqrt(1^2 + 1^2), sqrt(2^2 + 1^2)
%                         inf, sqrt(1^2 + 2^2), sqrt(0^2 + 2^2), sqrt(1^2 + 2^2), inf
%                        ];
      
% 3x3 neighbor grid, using Euclidean distance as cost            
% neighbor_cost_matrix = [sqrt(1^2 + 1^2), 1, sqrt(1^2 + 1^2)
%                         1, inf, 1
%                         sqrt(1^2 + 1^2), 1, sqrt(1^2 + 1^2)                      
%                        ];                   
             
[row,col] = ind2sub(size(neighbor_cost_matrix), find(~isinf(neighbor_cost_matrix)));     
transition_cost = neighbor_cost_matrix(find(~isinf(neighbor_cost_matrix)));
kernel_center = (((size(neighbor_cost_matrix)-[1,1])/2)+[1,1]);
neighbors = [col - kernel_center(2), ... % x
             row - kernel_center(1)];    % y

            
            
%% 
cost = inf * ones(size(track,1), size(track,2));
came_from_position = {};
came_from_position{size(track,1), size(track,2)} = [];
visited = zeros(size(track,1), size(track,2));

% Starting index
start = [160, 200];
goal = [135, 230];
figure(1);
immarker(start(1), start(2), 4, 'rx', 1);
immarker(goal(1), goal(2), 4, 'gx', 1);

% Set cost to 0 at starting node
cost(start(2), start(1)) = 0;

% Create priority queue and add starting point to it
q = PriorityQueue();
%q.push(start, 0); % starting cost/priority set to 0
q.insert([0, start]); % starting cost/priority set to 0

% ind2sub
% sub2ind(size(track), 1, 2)
iterations = 0;

GoalReached = false;
tic
%while (~q.empty())
while (q.size() > 0)
    % Pop next prioritized position to process
    %current = q.pop(); 
    current = q.remove();
    current = current(2:3);
    current_cost = cost(current(2), current(1));
    
    % Have we reached the goal?
    if isequal(current, goal)     
        GoalReached = true;
        disp('Goal reached');
        break;
    end
    
    % Evaluate all neighbors
    for i = 1:size(neighbors,1)
        neighbor = current + neighbors(i, :);
        
        % Verify that the neighbor is valid
        if (neighbor(1) > 0 && neighbor(1) <= size(track,2) && ...
            neighbor(2) > 0 && neighbor(2) <= size(track,1) && ...
            track(neighbor(2), neighbor(1)) == 1)
            visited(neighbor(2), neighbor(1)) = 1;
            
            dist_to_goal = sqrt(sum((goal - neighbor) .^ 2));            
        
            % Compute new cost for going to neigbor
            new_cost = current_cost + transition_cost(i);
            % If the new cost is better than the previously achieved cost for
            % the neighbor position, then update it
            if new_cost < cost(neighbor(2), neighbor(1))
                % Update cost
                cost(neighbor(2), neighbor(1)) = new_cost;
                came_from_position{neighbor(2), neighbor(1)} = current;
                
                % Add neighbor position to prioritized processing
                priority = new_cost + 0*dist_to_goal; % adding this dist_to_goal makes it an A* algorithm                
                %q.push(neighbor, priority); % could also be push_update (which checks if the element already exists)
                q.insert([priority, neighbor]); % could also be push_update (which checks if the element already exists)
            end
        end
    end
        
    iterations = iterations + 1;
    
    if mod(iterations, 100) == 0
        R = 255 * track;
        G = 255 * xor(track, visited);
        B = 255 * xor(track, visited);
        im(:,:,1) = R;
        im(:,:,2) = G;
        im(:,:,3) = B;
        imshow(im);    
        pause(0.001);
    end
end
toc

R = 255 * track;
G = 255 * xor(track, visited);
B = 255 * xor(track, visited);
im(:,:,1) = R;
im(:,:,2) = G;
im(:,:,3) = B;
imshow(im);    
pause(0.001);

if (~GoalReached)
    error('The goal was not reached');
end

% Backtrace path by starting in goal location
current = goal;
path = [current];
impath = zeros(size(track));
impath(current(2), current(1)) = 1;
while ~isequal(current, start)
    current = came_from_position{current(2), current(1)};
    path(end+1, :) = current;
    impath(current(2), current(1)) = 1;
end
path = path(end:-1:1,:);


R = 255 * track;
G = 255 * xor(track, impath);
B = 255 * xor(track, impath);
im(:,:,1) = R;
im(:,:,2) = G;
im(:,:,3) = B;
figure(2);
imshow(im);   

%%
guidance = 0 * ones(size(track,1), size(track,2));
guidance(sub2ind(size(track), path(:,1), path(:,2))) = 1;
guidance = -double(bwdist(1-track,'euclidean'));
imshow(1-guidance / min(min(guidance)));