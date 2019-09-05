function corners = getBoundingBox(pcd)
    min_x = inf;
    max_x = -inf;
    min_y = inf;
    max_y = -inf;
    min_z = inf;
    max_z = -inf;
        
    for (i = 1:length(pcd.x))
        if (pcd.x(i) > max_x)
            max_x = pcd.x(i);
        end
        if (pcd.x(i) < min_x)
            min_x = pcd.x(i);
        end
        if (pcd.y(i) > max_y)
            max_y = pcd.y(i);
        end
        if (pcd.y(i) < min_y)
            min_y = pcd.y(i);
        end        
        if (pcd.z(i) > max_z)
            max_z = pcd.z(i);
        end
        if (pcd.z(i) < min_z)
            min_z = pcd.z(i);
        end            
    end
    
    corners = [min_x, min_y, min_z;
               min_x, min_y, max_z;
               min_x, max_y, min_z;
               min_x, max_y, max_z;
               max_x, min_y, min_z;
               max_x, min_y, max_z;
               max_x, max_y, min_z;
               max_x, max_y, max_z];
end