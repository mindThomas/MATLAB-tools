function newPCD = removeFromPCD(pcd, idx)    
    newPCD = pcd;
    
    newPCD.x(idx) = [];
    newPCD.y(idx) = [];
    newPCD.z(idx) = [];
end