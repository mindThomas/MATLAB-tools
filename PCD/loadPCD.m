function pcd = loadPCD(filepath)
    tmp = loadpcd(filepath);

    % parse the PCD
    % x y z
    pcd.x = tmp(1,:)';
    pcd.y = tmp(2,:)';
    pcd.z = tmp(3,:)';