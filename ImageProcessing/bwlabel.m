function [ labels ] = bwlabel( data, varargin )
% bwlabel binary image labeling
% Labels a binary image through 8-point connectivity without the need for any toolboxes. 

if (nargin == 2)    
    connectivity = varargin{1};
    if (connectivity ~= 4 && connectivity ~= 8)
        error('Incorrect connectivity. Can only be 4 or 8');
    end
else
    connectivity = 8;
end

[x,y] = size(data);
% expand dataset to avoid crash when searching:
data = [zeros(1,y+2);[zeros(x,1) data zeros(x,1)]];
[x,y] = size(data);

labels = zeros(size(data));
nextlabel = 1;
linked = [];

for i = 2:x                       % for each row
    for j = 2:y-1                 % for each column
        if data(i,j) ~= 0         % not background
            % find binary value of neighbours
            if (connectivity == 8)
                neighboursearch = [data(i-1,j-1), data(i-1,j), data(i-1,j+1),data(i,j-1)];
            elseif (connectivity == 4)
                neighboursearch = [data(i-1,j),data(i,j-1)];
            end
            
            % search for neighbours with binary value 1
            [~,n,neighbours] = find(neighboursearch==1);
            
            % if no neighbour is allready labeled: assign new label
            if isempty(neighbours)
                linked{nextlabel} = nextlabel;
                labels(i,j) = nextlabel;
                nextlabel = nextlabel+1;                
            
            % if neighbours is labeled: pick the lowest label and store the
            % connected labels in "linked"
            else
                if (connectivity == 8)
                    neighboursearch_label = [labels(i-1,j-1), labels(i-1,j), labels(i-1,j+1),labels(i,j-1)];
                elseif (connectivity == 4)
                    neighboursearch_label = [labels(i-1,j), labels(i,j-1)];
                end
                L = neighboursearch_label(n);
                labels(i,j) = min(L);
                for k = 1:length(L)
                    label = L(k);
                    linked{label} = unique([linked{label} L]);
                end                
            end
        end
    end
end

% remove the previous expansion of the image
labels = labels(2:end,2:end-1);


%% join linked areas
% for each link, look through the other links and look for common labels.
% if common labels exist they are linked -> replace both link with the 
% union of the two. Repeat until there is no change in the links.

change2 = 1;
while change2 == 1
    change = 0;
    for i = 1:length(linked)
        for j = 1:length(linked)
            if i ~= j
                if sum(ismember(linked{i},linked{j}))>0 && sum(ismember(linked{i},linked{j})) ~= length(linked{i})
                    change = 1;
                    linked{i} = union(linked{i},linked{j});
                    linked{j} = linked{i};
                end
            end
        end
    end
    
    if change == 0
        change2 = 0;
    end
    
end

% removing redundant links
linked = unique(cellfun(@num2str,linked,'UniformOutput',false));
linked = cellfun(@str2num,linked,'UniformOutput',false);

K = length(linked);
templabels = labels;
labels = zeros(size(labels));

% label linked labels with a single label
for k = 1:K
    for l = 1:length(linked{k})
        labels(templabels == linked{k}(l)) = k;
    end
end

end

            
            
            
