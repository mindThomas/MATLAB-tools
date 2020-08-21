function pcd = loadPCD(filepath, varargin)
    pcd = loadpcd_auto(filepath);
    
    if (length(varargin) == 1)
        minimalLoad = varargin{1};
    else
        minimalLoad = false;
    end    

    if (minimalLoad == false)
        pcd.sensor = pcdHeader(filepath);
        pcd.sensor.tf = [pcd.sensor.rotm, pcd.sensor.origin; 0,0,0,1];
    end
