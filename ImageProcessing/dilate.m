function out = dilate(image)
    % structuring element
    s = [1 1 1;
         1 1 1;
         1 1 1];

    in = false([size(image,1)+size(s,1)-1, size(image,2)+size(s,2)-1]);
    in(2:end-1, 2:end-1) = image;
    out_padded = false(size(in));

    for (x = 2:(size(in,2)-1))
        for (y = 2:(size(in,1)-1))
            out_padded(y,x) = any( s & in((y-1):(y+1),(x-1):(x+1) ), 'all'); % or operation
        end
    end

    out = out_padded(2:end-1, 2:end-1);
end