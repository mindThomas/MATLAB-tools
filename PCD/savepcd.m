function savepcd(pcd, fname)   
    
    fp = fopen(fname, 'w');
        
    % unorganized point cloud
    npoints = length(pcd.x);
    width = npoints;
    height  = 1;    

    fields = 'x y z intensity';          
    siz = '4 4 4 1';
    typ = 'F F F U';  
    count = '1 1 1 1';    
    
    % write the PCD file header    
    fprintf(fp, '# .PCD v.7 - Point Cloud Data file format\n');
    fprintf(fp, 'VERSION .7\n');
    
    fprintf(fp, 'FIELDS %s\n', fields);
    fprintf(fp, 'SIZE %s\n', siz);
    fprintf(fp, 'TYPE %s\n', typ);
    fprintf(fp, 'COUNT %s\n', count);
        
    fprintf(fp, 'WIDTH %d\n', width);
    fprintf(fp, 'HEIGHT %d\n', height);
    fprintf(fp, 'POINTS %d\n', npoints);       
    
    % Pointcloud data
    data = [];
    data(1,:) = pcd.x;
    data(2,:) = pcd.y;
    data(3,:) = pcd.z;
    data(4,:) = pcd.intensity;
        
    % Write ASCII format data
    fprintf(fp, 'DATA ascii\n');
    fprintf(fp, '%f %f %f %d\n', data);    

    fclose(fp);
end


    

