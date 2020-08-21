function newPCD = extractPCD(pcd,idx)    
    newPCD.x = pcd.x(idx);
    newPCD.y = pcd.y(idx);
    newPCD.z = pcd.z(idx);         

    if (isfield(pcd, 'sensor'))
        newPCD.sensor = pcd.sensor;
    end    
    
    fields_in = fieldnames(pcd); 
    fields_newPCD = fieldnames(newPCD);      
    missingFields = setxor(fields_in, fields_newPCD);
    
    for (i = 1:length(missingFields))
        newPCD.(missingFields{i}) = pcd.(missingFields{i})(idx);
    end       
end
