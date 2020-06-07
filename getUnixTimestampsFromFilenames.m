function timestamps = getUnixTimestampsFromFilenames(files)
    timestamps = [];
    for (i = 1:length(files))
        filePath = files{i};
        %file = dir(filePath);
        [~, fName, ext] = fileparts(filePath);
        sName = split(fName, '_');
        fName = [sName{end-1}, '_', sName{end}];
        d = datetime(fName, 'InputFormat', 'yyyy-MM-dd_HH-mm-ss-SSS');
        timestamps(i) = posixtime(d);
    end
end
