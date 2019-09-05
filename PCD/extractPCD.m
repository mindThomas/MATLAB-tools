function newPCD = extractPCD(pcd,idx)    
    newPCD.x = pcd.x(idx);
    newPCD.y = pcd.y(idx);
    newPCD.z = pcd.z(idx);
    if (isfield(pcd, 'sensor'))
        newPCD.sensor = pcd.sensor;
    end
end