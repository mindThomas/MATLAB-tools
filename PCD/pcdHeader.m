function [sensor] = pcdHeader(fname)
    
    fp = fopen(fname, 'r');
    version = [];

    while true
        line = fgetl(fp);
        
        if line(1) == '#'
            continue;
        end
        
        [field,remain] = strtok(line, ' \t');
        remain = strtrim(remain);
        
        switch field
            case 'VERSION'
                version = remain;
            case 'FIELDS'
                fields = remain;
            case 'TYPE'
                type = remain;
            case 'WIDTH'
                width = str2num(remain);
            case 'HEIGHT'
                height = str2num(remain);
            case 'POINTS'
                npoints = str2num(remain);
            case 'SIZE'
                siz = str2num(remain);
            case 'COUNT'
                count = str2num(remain);
            case 'VIEWPOINT'
                viewpoint = str2num(remain);                
            case 'DATA'
                mode = remain;
                break;
            otherwise
                fprintf('unknown field %s\n', field);
        end
    end
    
    % if no version field we'll assume it's not a PCD file
    if isempty(version)
        return;
    end
    
    if height == 1
        org = 'unorganized';
    else
        org = 'organized';
    end
    %fprintf('%s: %s, %s, <%s> %dx%d\n', ...
    %    fname, mode, org, fields, width, height);
    %fprintf('  %s; %s\n', type, num2str(siz));
    
    sensor.origin = viewpoint(1:3)';
    sensor.quaternion = viewpoint(4:7)';    
    sensor.rotm = quat2rotm(sensor.quaternion);
    
    fclose(fp);
end