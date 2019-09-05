function out = rotatePCD(pcd, R)
    points = [pcd.x, pcd.y, pcd.z];
    rotated = (R * points')';
    out = pcd;
    out.x = rotated(:,1);
    out.y = rotated(:,2);
    out.z = rotated(:,3);
end