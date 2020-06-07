function [timestamp, unix_seconds] = timestampFromFilename(file)    
    [filedir,s,ext] = fileparts(file);

    timestamp = 0;
    underscores = find(s == '_');
    if (length(underscores) >= 2)
        s2 = s(1:underscores(2)-1);
    else    
        dot = find(s == '.');
        if (length(dot) == 1)
            s2 = s(1:dot(1)-1);
        else
            return;
        end
    end
    
    timestamp = datetime(s2, 'InputFormat', 'yyyy-MM-dd_HH-mm-ss-SSS');
    unix_seconds = posixtime(timestamp);
end
