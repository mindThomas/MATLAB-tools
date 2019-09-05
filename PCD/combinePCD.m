function out = combinePCD(pcd1, pcd2)
    % Compute transform between the two pointclouds
    pcd1_to_world = [pcd1.sensor.rotm, pcd1.sensor.origin; 0,0,0,1];
    pcd2_to_world = [pcd2.sensor.rotm, pcd2.sensor.origin; 0,0,0,1];   
    pcd2_to_pcd1 = inv(pcd1_to_world) * pcd2_to_world;
    
    % Transform pcd2 points into pcd1 frame
    pcd2_transformed = transformPCD(pcd2, pcd2_to_pcd1);
    
    % Store combined output    
    out.x = [pcd1.x; pcd2_transformed.x];
    out.y = [pcd1.y; pcd2_transformed.y];
    out.z = [pcd1.z; pcd2_transformed.z];
    
    % Set output transform to pcd1
    out.sensor = pcd1.sensor; 
end