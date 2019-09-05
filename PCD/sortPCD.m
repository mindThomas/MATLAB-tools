function sorted = sortPCD(pcd)
    points = [pcd.x,pcd.y];
    distance = sqrt(sum(points .* points, 2));
    [sortedDist, sortedIdx] = sort(distance, 'ascend');
    
    sorted = extractPCD(pcd, sortedIdx);
end