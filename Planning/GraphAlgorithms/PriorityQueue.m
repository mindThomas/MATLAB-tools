classdef PriorityQueue < matlab.mixin.Copyable % could also have been "handle"
    % https://www.redblobgames.com/pathfinding/a-star/introduction.html
	properties %(SetAccess = private)
        
    end  
    properties (SetAccess = private)
        elements % elements such as index
        priority % priority weight (lower value = higher priority)
    end
    methods
        function obj = PriorityQueue()
            obj.elements = {};
            obj.priority = []; 
        end        
        
        function push(obj, element, priority)
            % Push back a priority element
            obj.elements{end+1} = element;
            obj.priority(end+1) = priority;
        end
        
        function push_update(obj, element, priority)            
            % Push back a priority element only if it doesn't exist in the
            % queue already
            %element_indices = find(obj.elements == element);
            element_indices = find(cellfun(@(x) isequal(x, element), obj.elements));
            if ~isempty(element_indices)
                % Element already exists, so update its priority                
                obj.priority(element_indices) = priority;    
            else
                % Add new element                
                obj.elements{end+1} = element;
                obj.priority(end+1) = priority;
            end
        end

        function element = pop(obj)
            % Pop element with highest priority (lower weight)
            [~, queue_index] = min(obj.priority);
            if ~isempty(queue_index)
                element = obj.elements{queue_index}; % get stored element index
            
                % Remove element from queue
                obj.elements(queue_index) = [];
                obj.priority(queue_index) = [];
            else
                element = {};
            end
        end
        
        function update_priority(obj, element, priority)
            %element_indices = find(obj.elements == element);
            element_indices = find(cellfun(@(x) isequal(x, element), obj.elements));
            if ~isempty(element_indices)
                obj.priority(element_indices) = priority;            
            end
        end      
        
        function is_empty = empty(obj)
           is_empty = isempty(obj.priority);
        end
    end
end