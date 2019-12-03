classdef OccupancyGrid
    properties
        pmin % upper left corner position [x,y]
        pmax % lower right corner position [x,y]
        resolution % x-y cell size
        grid % probability of cell being occupied
        p0
        l0   % initial log odds (initial probability of grid cells)
    end
    
    methods
        function obj = OccupancyGrid(pmin, pmax, resolution)
            obj.pmin = pmin;
            obj.pmax = pmax;
            obj.resolution = resolution;
                        
            cell_count = ceil((pmax-pmin) ./ resolution);
            obj.resolution = (pmax-pmin) ./ cell_count;            
            cell_count = cell_count + ones(size(pmin));
            
            obj.p0 = 0.5;
            obj.grid = obj.p0 * ones(cell_count(2), cell_count(1)); % set all cells to probability of 0.5 of being occupied
            obj.l0 = obj.logOdds(obj.p0); % log odds of initial probability
        end
        
        function obj = load(obj, filename)
            im = imread(filename);
            if (size(im, 1) ~= size(obj.grid, 1) || size(im, 2) ~= size(obj.grid, 2))
                error('Incorrect map resolution');
            end
            obj.grid = 1-im;
        end
        
        function obj = RangeBearing_Update(obj, likelihood_field)
            % Using a likelihood field computed from the inverse 
            % measurement model using the range-bearing measurement model
            % (see BeamRangeSensor.m) based on Probabilistic Robotics
            %
            % This function incorporates a full likelihood field map
            % (should be of the same size and with same origin as this
            % occupancy grid) into the current grid using Binary Bayes
            % filter for the update. To handle the extremeties the log odds
            % are used to handle the probabilities of each grid cell
            %
            % See Table 9.1 in Probabilistic Robotics
            %logOdds = obj.toLogOdds(obj.grid);
            
            % Loop over all cells and add likelihood
            for (i = 1:size(obj.grid, 1))
                for (j = 1:size(obj.grid, 2))                  
                    if (likelihood_field(i,j) ~= obj.p0)
                        obj.grid(i,j) = obj.probabilityFromLogOdds( obj.logOdds(obj.grid(i,j)) + obj.logOdds(likelihood_field(i,j)) - obj.l0 );
                    end
                    %logOdds(i,j) = logOdds(i,j) + obj.logOdds(likelihood_field(i,j)) - obj.l0;
                end
            end  
            
            % Convert back into probability map            
            %obj.grid = obj.fromLogOdds(logOdds);
        end
        
        function obj = Update_MaximumAposteriori(obj, scan, pose)
            % See Table 9.3 in Probabilistic Robotics
            % This implements MAP estimation of the map posterior
            % Thereby yielding the map that leads to the highest posterior
            % probability/likelihood
            % See especially figure 9.11 for a example
            % (compare 9.11b to 9.11d)
        end
        
        function s = size(obj)
            s = size(obj.grid);
        end
        
        function [x,y] = getXY(obj, row, column)
            x = obj.pmin(1) + obj.resolution(1)*(column-1);
            y = obj.pmax(2) - obj.resolution(2)*(row-1);
        end
        
        function [row,col] = getRowCol(obj, x, y)
            col = round((x - obj.pmin(1)) / obj.resolution(1)) + 1;
            row = round((obj.pmax(2) - y) / obj.resolution(2)) + 1;
        end
        
        function [x,y] = getXYfromIndex(obj, idx)
            [row,col] = ind2sub(size(obj.grid), idx);
            [x,y] = obj.getXY(row, col);
        end
        
        function plot(obj)
            % Visualize the occupancy grid in grayscale
            % Invert the probabilities such that free cells, which has
            % probability values of 0.0 (of being occupied), should be
            % colored as white
            im = ones(size(obj.grid)) - obj.grid;
            imshow(im);
            %xticks([0 10 20 30 40 50])
            axis on;
            tickres = (obj.pmax - obj.pmin) / 10;
            tickstep = round(tickres ./ obj.resolution);
            tickres = tickstep .* obj.resolution;
            
            t = 1:tickstep(1):size(obj.grid,2);
            xticks(t);
            
            t = 1:tickstep(2):size(obj.grid,1);
            yticks(t);
            
            t = (obj.pmin(1):tickres(1):obj.pmax(1))';
            x_ticks = t;
            
            t = (obj.pmin(2):tickres(2):obj.pmax(2))';
            y_ticks = t(end:-1:1);
            xticklabels(num2str(x_ticks));
            yticklabels(num2str(y_ticks));
            
            %xticklabel
            %xlim([obj.pmin(1), obj.pmax(1)]);
            %ylim([obj.pmin(2), obj.pmax(2)]);
        end
        
        function oddsMap = toOdds(obj, probability_map)
            % The odds of a state x is defined as the ratio of the
            % probability of this event divided by the probability of its negate
            % odds = p / (1 - p);
            
            oddsMap = zeros(size(probability_map));
            % Loop over all cells and compute odds
            for (i = 1:size(probability_map, 1))
                for (j = 1:size(probability_map, 2))
                    p = probability_map(i,j);
                    oddsMap(i,j) = p / (1 - p);
                end
            end            
        end   
        
        function lodds = logOdds(obj, p)
            if (p > 1)
                p = 1;
            end
            if (p < 0)
                p = 0;
            end
            lodds = log( p / (1 - p) );
        end
        
        function p = probabilityFromLogOdds(obj, logOdds)
            p = 1 - 1 / (1+exp(logOdds));
        end
        
        function logOddsMap = toLogOdds(obj, probability_map)
            % For binary states the belief (posterior of the state) is
            % commonly implemented as a log odds ratio.
            % Log odds assume values from -inf to inf
            % l = log( p / (1 - p) );
            
            logOddsMap = zeros(size(probability_map));
            % Loop over all cells and compute log odds
            for (i = 1:size(probability_map, 1))
                for (j = 1:size(probability_map, 2))
                    p = probability_map(i,j);
                    logOddsMap(i,j) = log( p / (1 - p) );
                end
            end              
        end 
        
        function probability_map = fromLogOdds(obj, logOddsMap)
            % Convert log odds to probability probability of the state being 1
            % p(x) = p(x=1)
            % p = 1 - 1 / (1+exp(obj.l));            
            
            probability_map = zeros(size(logOddsMap));
            % Loop over all cells and compute log odds
            for (i = 1:size(logOddsMap, 1))
                for (j = 1:size(logOddsMap, 2))
                    p = 1 - 1 / (1+exp(logOddsMap(i,j)));
                    probability_map(i,j) = p;
                end
            end             
        end        
    end
end