global ds track robot_length robot_width;
robot_length = 1.7;
robot_width = 1;
ds = 0.1; % grid resolution

goal_x = 9;
goal_y = 25;

raw_track = imread('track.png');
track = raw_track(:,:,1) > 0;

figure(1);
image([0 ds*size(track,2)],[0 ds*size(track,1)], repmat(flip(track), [1,1,3])); %, 'AlphaData', 1.0);
set(gca,'YDir','normal');
xlim([0 ds*size(track,2)]);
ylim([0 ds*size(track,1)]);
axis equal;

%% Create distance transform from goal position
% Create neighbour cost
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
transition_distance = neighbor_cost_matrix(find(~isinf(neighbor_cost_matrix)));
kernel_center = (((size(neighbor_cost_matrix)-[1,1])/2)+[1,1]);
neighbors = [col - kernel_center(2), ... % x
             row - kernel_center(1)];    % y

% Create the distance image
distance = inf * ones(size(track));
        
% Assuming bottom left corner of grid is (0,0)
goal_x_idx = round(goal_x / ds) + 1;
goal_y_idx = size(distance, 1) - round(goal_y / ds);

% Create priority queue and add starting point to it
q = PriorityQueue();
%q.push(start, 0); % starting cost/priority set to 0
q.insert([0, goal_x_idx, goal_y_idx]); % starting cost/priority set to 0

while (q.size() > 0)
    % Pop next prioritized position to process
    %current = q.pop(); 
    element = q.remove();
    
    current_distance = element(1);
    current = element(2:3);    
    
    % Evaluate all neighbors
    for i = 1:size(neighbors,1)
        neighbor = current + neighbors(i, :);
        
        % Verify that the neighbor is valid
        if (neighbor(1) > 0 && neighbor(1) <= size(distance,2) && ...
            neighbor(2) > 0 && neighbor(2) <= size(distance,1) && ...
            track(neighbor(2), neighbor(1)) == 1)
                    
            % Compute new cost for going to neigbor
            new_distance = current_distance + transition_distance(i);
            % If the new cost is better than the previously achieved cost for
            % the neighbor position, then update it
            if new_distance < distance(neighbor(2), neighbor(1))
                % Update cost
                distance(neighbor(2), neighbor(1)) = new_distance;                
                
                % Add neighbor position to prioritized processing
                %q.push(neighbor, priority); % could also be push_update (which checks if the element already exists)
                q.insert([new_distance, neighbor]); % could also be push_update (which checks if the element already exists)
            end
        end
    end       
end

guidance_distance = distance;

%%
% Robot init pose
p = [18, 34, deg2rad(45)]; % [x, y, heading]

r = Robot(robot_length, robot_width, p);
hold on;
r.draw();
hold off;

%%
% Create priority queue
% Queue element = [exploration_cost (priority), cost, x, y, theta]
q = PriorityQueue(1);
q.insert([0, 0, p]);

curvatures = -0.4:0.1:0.4;
distances = 0.5:0.1:1;
while (q.size() > 0)
    element = q.remove();
        
    p = element(3:5);
    current_cost = p(2);

    [A,B] = meshgrid(curvatures, distances);
    actuation_guess = [A(:), B(:)]; % [curvature, distance]

    poses = Robot.step(p, actuation_guess(:,1), actuation_guess(:,2));

    figure(1);
    hold on;
    Robot.draw_(p(1), p(2), p(3), robot_length, robot_width, 'r');
    hold off;
    
    redraw(p);
    hold on;
    proposals = 0;
    for (i = 1:size(poses, 1))
        collision = Robot.in_collision(poses(i,1), poses(i,2), poses(i,3), r.length, r.width, track, ds);
        if (~collision)
            Robot.draw_(poses(i,1), poses(i,2), poses(i,3), r.length, r.width);
            
            % Assuming bottom left corner of grid is (0,0)
            ix = round(poses(i,1) / ds) + 1;
            iy = size(guidance_dist, 1) - round(poses(i,2) / ds);
            
            % Add proposal
            distance = actuation_guess(i, 2);            
            new_cost = current_cost + distance;
            exploration_cost = new_cost + 10*guidance_distance(iy,ix);
            q.insert([exploration_cost, new_cost, poses(i,:)]);
        end    
    end
    hold off;
    drawnow;
end


function redraw(p)
    global ds track robot_length robot_width;
    figure(2);
    image([0 ds*size(track,2)],[0 ds*size(track,1)], repmat(flip(track), [1,1,3])); %, 'AlphaData', 1.0);
    set(gca,'YDir','normal');
    xlim([0 ds*size(track,2)]);
    ylim([0 ds*size(track,1)]);
    axis equal;
    
    hold on;
    Robot.draw_(p(1), p(2), p(3), robot_length, robot_width, 'r');
    hold off;
end