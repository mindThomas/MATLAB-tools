classdef Line2
   properties
      a
      b
      c
   end
   properties (Access = private)
       ShouldIntersectOrigin
   end
   methods
      function obj = Line2(ShouldIntersectOrigin_, a, b, c)
         if (nargin == 1)
            obj.ShouldIntersectOrigin = ShouldIntersectOrigin_;
            obj.a = 0;
            obj.b = 0;            
            obj.c = 0;   
         elseif (nargin == 4)
            obj.ShouldIntersectOrigin = ShouldIntersectOrigin_;
            obj.a = a;
            obj.b = b;            
            obj.c = c;            
         else
            obj.ShouldIntersectOrigin = false;
            obj.a = 0;
            obj.b = 0;            
            obj.c = 0;   
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
         z = obj.a*x + obj.b*x.^2 + obj.c;
      end
      
      function plot(obj, x_min, x_max)          
          fplot(@(x) obj.eval(x), [x_min, x_max]);
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
            A = [x, x.^2];
         else
            A = [x, x.^2, ones(length(x),1)];
         end
         X = lscov(A, Y);

         obj.a = X(1);
         obj.b = X(2);
         if (obj.ShouldIntersectOrigin)
            obj.c = 0;
         else
            obj.c = X(3);
         end
      end

      % Weighted Least Squares fitting of a regular 2D plane
      function obj = fitWLS(obj, x, y, w)
         % y = A*x + b
         % Putting b into the parameter vector, x, by adding a column of 1's to A
         % y = A*x
         Y = y;
         if (obj.ShouldIntersectOrigin)
            A = [x, x.^2];
         else
            A = [x, x.^2, ones(length(x),1)];
         end
         X = lscov(A, Y, w);

         obj.a = X(1);
         obj.b = X(2);
         if (obj.ShouldIntersectOrigin)
            obj.c = 0;
         else
            obj.c = X(3);
         end
      end
   end
end
