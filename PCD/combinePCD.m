function out = combinePCD(pcd1, pcd2, varargin)    
    if (length(varargin) == 1)
        enableTransform = varargin{1};
    else
        enableTransform = true;
    end
    % Compute transform between the two pointclouds
    pcd1_to_world = [pcd1.sensor.rotm, pcd1.sensor.origin; 0,0,0,1];
    pcd2_to_world = [pcd2.sensor.rotm, pcd2.sensor.origin; 0,0,0,1];   
    pcd2_to_pcd1 = inv(pcd1_to_world) * pcd2_to_world;
    
    if (enableTransform)
        % Transform pcd2 points into pcd1 frame
        pcd2_transformed = transformPCD(pcd2, pcd2_to_pcd1);
    else
        pcd2_transformed = pcd2;
    end
    
    % Store combined output    
    if (isfield(pcd1, 'x'))
        out.x = [pcd1.x; pcd2_transformed.x];
    end
    if (isfield(pcd1, 'y'))
        out.y = [pcd1.y; pcd2_transformed.y];
    end
    if (isfield(pcd1, 'z'))
        out.z = [pcd1.z; pcd2_transformed.z];
    end      

    % Set output transform to pcd1
    out.sensor = pcd1.sensor;   
    
    fields_in = fieldnames(pcd1); 
    fields_out = fieldnames(out);      
    missingFields = setxor(fields_in, fields_out);
    
    for (i = 1:length(missingFields))
        if (isfield(pcd1, missingFields{i}) && isfield(pcd2, missingFields{i}))
            out.(missingFields{i}) = [pcd1.(missingFields{i}); pcd2.(missingFields{i})]; 
        end
    end             
end
