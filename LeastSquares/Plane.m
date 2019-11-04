classdef Plane
   properties
      a
      b
      c
      sigma2
   end
   methods
      function obj = Plane(a, b, c)
         if (nargin == 3)
            obj.a = a;
            obj.b = b;
            obj.c = c;
         else
            obj.a = 0;
            obj.b = 0;
            obj.c = 0;
         end
      end

      function obj = fit(obj, x, y, z, varargin)
         if (length(varargin) == 1)
            obj = fitWLS(obj, x, y, z, varargin{1});
         else
            obj = fitOLS(obj, x, y, z);
         end
      end

      function z = eval(obj, x, y)
         z = obj.a*x + obj.b*y + obj.c;
      end
   end

   methods (Access = private)     
      % Ordinary Least Squares fitting of a regular 2D plane
      function obj = fitOLS(obj, x, y, z)
         % y = A*x + b
         % Putting b into the parameter vector, x, by adding a column of 1's to A
         % y = A*x
         Y = z;
         A = [x, y, ones(length(x),1)];
         X = lscov(A, Y);

         obj.a = X(1);
         obj.b = X(2);
         obj.c = X(3);
      end

      % Weighted Least Squares fitting of a regular 2D plane
      function obj = fitWLS(obj, x, y, z, w)
         % y = A*x + b
         % Putting b into the parameter vector, x, by adding a column of 1's to A
         % y = A*x
         Y = z;
         A = [x, y, ones(length(x),1)];
         X = lscov(A, Y, w);

         obj.a = X(1);
         obj.b = X(2);
         obj.c = X(3);
      end
   end
end
