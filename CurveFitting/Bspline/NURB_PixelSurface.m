classdef NURB_PixelSurface
    % Takes in a greyscale image and sets up a weighted 3-order 2D basis function on
    % each pixel using the greyscale value as the weight
    properties %(SetAccess = private)
        image
    end
    methods
        function obj = NURB_PixelSurface(image)                                    
            % pad image with two rows/columns of pixels around all edges
            obj.image = zeros(size(image)+[4,4]);
            obj.image(3:end-2, 3:end-2) = double(image);
            
            % Duplicate pixels near the edges such that the value doesn't
            % go to 0
            obj.image(3:end-2, 1) = obj.image(3:end-2, 3);
            obj.image(3:end-2, 2) = obj.image(3:end-2, 3);            
            
            obj.image(3:end-2, end) = obj.image(3:end-2, end-2);
            obj.image(3:end-2, end-1) = obj.image(3:end-2, end-2);
            
            obj.image(1, 3:end-2) = obj.image(3, 3:end-2);
            obj.image(2, 3:end-2) = obj.image(3, 3:end-2);
            
            obj.image(end, 3:end-2) = obj.image(end-2, 3:end-2);
            obj.image(end-1, 3:end-2) = obj.image(end-2, 3:end-2);
            
            % Duplicate the corner pixels
            obj.image(1:2, 1:2) = obj.image(3, 3);                        
            obj.image(end-1:end, 1:2) = obj.image(end-2, 3);                        
            obj.image(end-1:end, end-1:end) = obj.image(end-2, end-2);                        
            obj.image(1:2, end-1:end) = obj.image(3, end-2);            
        end        
        
        function out = eval(obj, x, y)
            % Due to the choice of 3rd order basis function for a query points (x,y)
            % we need to consider all the control points in a neighbourhood
            % of:
            % x-2 <= x_i < x+2
            % y-2 <= y_i < y+2
            
            % First we find the corresponding nearest center
            x_c = floor(x);
            y_c = floor(y);
            x_min = x_c - 1;
            y_min = y_c - 1;
            x_max = x_c + 2;
            y_max = y_c + 2;
            
%             % Next we extract that 4x4 region of interest from the image
%             ROI = obj.image(y_min:y_max, x_min:x_max);
%       
%             % Now we sum up all the contributions
%             wROI = zeros(size(ROI));
%             for (xi = x_min:x_max)
%                 for (yi = y_min:y_max)
%                     wROI(yi-y_min+1, xi-x_min+1) = obj.basis(xi, x) * obj.basis(yi, y);
%                 end
%             end            
            
            % Now we sum up all the contributions
%             out = 0;
%             normalizer = 0;
%             for (xi = x_min:x_max)
%                 for (yi = y_min:y_max)
%                     activation = obj.basis(xi, x) * obj.basis(yi, y); 
%                     out = out + activation * obj.image(yi+2, xi+2);
%                     normalizer = normalizer + activation;                    
%                 end
%             end
%             out = out / normalizer;       
            % Normalization is not needed since the basis functions integrates and sums to 1  
%             out = 0;            
%             for (xi = x_min:x_max)
%                 for (yi = y_min:y_max)                    
%                     out = out + obj.image(yi+2, xi+2) * obj.basis(xi, x) * obj.basis(yi, y);                    
%                 end
%             end                            

            % Another way of doing this which leads to less evaluations 
            % of the basis functions is:
            out = 0;
            for (yi = y_min:y_max)    
                tmp = 0;        
                % tmp = get_row_value
                for (xi = x_min:x_max)
                    tmp = tmp + obj.image(yi+2, xi+2) * obj.basis(xi, x);                    
                end
                
                out = out + obj.basis(yi, y) * tmp;
            end                            

        end
        
        function B = basis(obj, i, u)       
            % Single centered 3-order basis function with uniformly distributed
            % knots, distributed with distance 1.
            % Basis function is centered around u=i
                         
            u_local = u - i;
            
            if (u_local >= -2 && u_local < -1)
                B = 1/6 * (u_local + 2)^3;
            elseif (u_local >= -1 && u_local < 0)
                B = -1/2*u_local^3 - u_local^2 + 2/3;
            elseif (u_local >= 0 && u_local < 1)
                B = 1/2*u_local^3 - u_local^2 + 2/3;
            elseif (u_local >= 1 && u_local < 2)
                B = 1/6 * (2 - u_local)^3;
            else
                B = 0;
            end
        end   
        
        function dB = dbasis(obj, i, u)       
            % Derivative of the 3-order basis above
            % Basis function is centered around u=i
                         
            u_local = u - i;
            
            if (u_local >= -2 && u_local < -1)
                dB = 3/6 * (u_local + 2)^2;
            elseif (u_local >= -1 && u_local < 0)
                dB = -3/2*u_local^2 - 2*u_local;
            elseif (u_local >= 0 && u_local < 1)
                dB = 3/2*u_local^2 - 2*u_local;
            elseif (u_local >= 1 && u_local < 2)
                dB = -3/6 * (2 - u_local)^2;
            else
                dB = 0;
            end
        end     
        
        function ddB = ddbasis(obj, i, u)       
            % Derivative of the 3-order basis above
            % Basis function is centered around u=i
                         
            u_local = u - i;
            
            if (u_local >= -2 && u_local < -1)
                ddB = (u_local + 2);
            elseif (u_local >= -1 && u_local < 0)
                ddB = -3*u_local - 2;
            elseif (u_local >= 0 && u_local < 1)
                ddB = 3*u_local - 2;
            elseif (u_local >= 1 && u_local < 2)
                ddB = (2 - u_local);
            else
                ddB = 0;
            end
        end       
   
        function dddB = dddbasis(obj, i, u)       
            % Derivative of the 3-order basis above
            % Basis function is centered around u=i
                         
            u_local = u - i;
            
            if (u_local >= -2 && u_local < -1)
                dddB = 1;
            elseif (u_local >= -1 && u_local < 0)
                dddB = -3;
            elseif (u_local >= 0 && u_local < 1)
                dddB = 3;
            elseif (u_local >= 1 && u_local < 2)
                dddB = -1;
            else
                dddB = 0;
            end
        end    
        
    end
end