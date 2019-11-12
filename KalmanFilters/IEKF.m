classdef EKF
	properties %(SetAccess = private)
        x        
    end    
	properties (SetAccess = private)        
        P
        K
        S
    end        
    properties (SetAccess = private)%(Access = private)
        f        
        h
        Q
        R           
        Q_sqrt
        R_sqrt 
        ts
    end
    properties (SetAccess = private)%(Access = private)
        JacobiansAvailable
        Fx
        Fu
        Fq
        Hx
        Hr
    end    
    methods
        function obj = IEKF(varargin)                        

        end        
     
        %% Base this on: https://iopscience.iop.org/article/10.1088/1742-6596/659/1/012022/pdf
end