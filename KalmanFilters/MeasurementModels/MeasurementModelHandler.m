classdef MeasurementModelHandler
	properties (Access = private)
        h
        R
        R_sqrt
    end
    methods
        function obj = MeasurementModelHandler(h, R)                                    
            obj.h = h;
            obj.R = R;
            obj.R_sqrt = chol(R, 'lower');
        end
        
        function z = get(obj, x)            
            % Generate a realization from the measurement model based on x and R
            % Generate noise realization, r[k] ~ N(0, R)
            r = obj.R_sqrt * randn([size(obj.R, 1), 1]);
            % Calculate the measurement model output
            z = obj.h(x, r);
        end
        
        function z = getDeterministic(obj, x)    
            % Generate a deterministic realization of the measurement model based on x and R
            % Zero noise
            r = zeros([size(obj.R, 1), 1]);
            % Calculate the measurement model output
            z = obj.h(x, r);
        end        
        
    end
end