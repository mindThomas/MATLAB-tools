classdef RunningAverage
    properties
        value
    end  
    properties (Access = private)
        n
    end
    
    methods
        function obj = RunningAverage()            
            obj.value = 0;
            obj.n = 0;
        end   
        
        function obj = Update(obj, in)
            prev_sum = obj.value * obj.n;
            new_sum = prev_sum + in;
            obj.n = obj.n + 1;
            obj.value = new_sum / obj.n;
        end
    end
end
