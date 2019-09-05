function newPCD = boxFilter(pcd, lim_min, lim_max)
    idx = find( (pcd.x >= lim_min(1)) & ...
                (pcd.x <= lim_max(1)) & ...
                (pcd.y >= lim_min(2)) & ...
                (pcd.y <= lim_max(2)) & ...
                (pcd.z >= lim_min(3)) & ...
                (pcd.z <= lim_max(3)) );
    
    newPCD = extractPCD(pcd, idx);
end