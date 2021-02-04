function out = Downsample(in, out_dt)

    t0 = in.time(1);
    indices = 1;
    for (i = 1:length(in.time))
        t1 = in.time(i);
        
        if ((t1-t0) > out_dt)
            t0 = t1;
            indices(end+1,1) = i;
        end
    end
    
    out = [];
    f = fieldnames(in);
    for i = 1:length(f)
        out.(f{i}) = in.(f{i})(indices,:);    
    end
    
end