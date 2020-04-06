classdef Line
   properties
      a
      b                  
   end
   properties (Access = private)
       ShouldIntersectOrigin
   end
   methods
      function obj = Line(ShouldIntersectOrigin_, a, b)
         if (nargin == 1)
            obj.ShouldIntersectOrigin = ShouldIntersectOrigin_;
            obj.a = 0;
            obj.b = 0;            
         elseif (nargin == 3)
            obj.ShouldIntersectOrigin = ShouldIntersectOrigin_;
            obj.a = a;
            obj.b = b;            
         else
            obj.ShouldIntersectOrigin = false;
            obj.a = 0;
            obj.b = 0;            
         end
      end

      function obj = fit(obj, x, y, varargin)
         if (length(varargin) == 1)
            obj = fitWLS(obj, x, y, varargin{1});
         else
            obj = fitOLS(obj, x, y);
         end
      end

      function z = eval(obj, x)
         z = obj.a*x + obj.b;
      end
      
      function plot(obj, x_min, x_max, varargin)
          if (length(varargin) == 0)
            line([x_min, x_max], [obj.eval(x_min), obj.eval(x_max)]);
          else
              line([x_min, x_max], [obj.eval(x_min), obj.eval(x_max)], 'Color', varargin{1});
          end
      end
   end

   methods (Access = private)     
      % Ordinary Least Squares fitting of a regular 2D plane
      function obj = fitOLS(obj, x, y)
         % y = A*x + b
         % Putting b into the parameter vector, x, by adding a column of 1's to A
         % y = A*x
         Y = y;
         if (obj.ShouldIntersectOrigin)
            A = [x];
         else
            A = [x, ones(length(x),1)];
         end
         X = lscov(A, Y);

         obj.a = X(1);
         if (obj.ShouldIntersectOrigin)
            obj.b = 0;
         else
            obj.b = X(2);         
         end
      end

      % Weighted Least Squares fitting of a regular 2D plane
      function obj = fitWLS(obj, x, y, w)
         % y = A*x + b
         % Putting b into the parameter vector, x, by adding a column of 1's to A
         % y = A*x
         Y = y;
         if (obj.ShouldIntersectOrigin)
            A = [x];
         else
            A = [x, ones(length(x),1)];
         end
         X = lscov(A, Y, w);

         obj.a = X(1);         
         if (obj.ShouldIntersectOrigin)
            obj.b = 0;
         else
            obj.b = X(2);         
         end
      end            
   end
end
