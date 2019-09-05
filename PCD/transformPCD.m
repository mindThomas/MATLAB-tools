function out = transformPCD(pcd, tf)
    points = [pcd.x, pcd.y, pcd.z, ones(length(pcd.x),1)];
    transformed = (tf * points')';
    out = pcd;
    out.x = transformed(:,1);
    out.y = transformed(:,2);
    out.z = transformed(:,3);
    
    if (isfield(pcd, 'sensor'))
        newOrigin = tf * [pcd.sensor.origin; 1;]
        out.sensor.origin = newOrigin(1:3);
        out.sensor.rotm = tf(1:3,1:3) * pcd.sensor.rotm;
    end
end