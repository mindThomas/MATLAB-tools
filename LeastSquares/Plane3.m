classdef Plane3
   properties
      a
      b
      c
      a2
      b2
      a3
      b3      
   end
   methods
      function obj = Plane3(a, b, c)
         if (nargin == 3)
            obj.a = a;
            obj.b = b;
            obj.c = c;
            obj.a2 = 0;
            obj.b2 = 0;
            obj.a3 = 0;
            obj.b3 = 0;
         else
            obj.a = 0;
            obj.b = 0;
            obj.c = 0;
            obj.a2 = 0;
            obj.b2 = 0;            
            obj.a3 = 0;
            obj.b3 = 0;            
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
         z = obj.a*x + obj.b*y + obj.a2*x.^2 + obj.b2*y.^2 + obj.a3*x.^3 + obj.b3*y.^3 + obj.c;
      end
   end

   methods (Access = private)    
      % Ordinary Least Squares fitting of a third-order 2D surface
      function obj = fitOLS(obj, x, y, z)
         % y = A*x + b
         % Putting b into the parameter vector, x, by adding a column of 1's to A
         % y = A*x
         Y = z;
         A = [x, y, x.^2, y.^2, x.^3, y.^3, ones(length(x),1)];
         X = lscov(A, Y);

         obj.a = X(1);
         obj.b = X(2);
         obj.a2 = X(3);
         obj.b2 = X(4);
         obj.a3 = X(5);
         obj.b3 = X(6);
         obj.c = X(7);         
      end

      % Weighted Least Squares fitting of a third-order 2D surface
      function obj = fitWLS(obj, x, y, z, w)
         % y = A*x + b
         % Putting b into the parameter vector, x, by adding a column of 1's to A
         % y = A*x
         Y = z;
         A = [x, y, x.^2, y.^2, x.^3, y.^3, ones(length(x),1)];

         % Downweigh large values of a2,a3,b2,b3 (resulting in large inclinations)
         n = length(x);         
         A = [A; 
              0,0,1,0,0,0,0;
              0,0,0,1,0,0,0;
              0,0,0,0,1,0,0;
              0,0,0,0,0,1,0];
         Y = [Y; 0;0;0;0];
         w = [w; 
              1000^2; % 0*n^2
              1000^2;
              1000^3; % 0*n^3
              1000^3];

         X = lscov(A, Y, w);

         obj.a = X(1);
         obj.b = X(2);
         obj.a2 = X(3);
         obj.b2 = X(4);
         obj.a3 = X(5);
         obj.b3 = X(6);
         obj.c = X(7); 
      end
   end
end
