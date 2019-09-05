function fex = fileExist(file_name)
    fex = false;

    dir_listing = dir(file_name);
    % check if "file_name" is a folder ...
    if (length(dir_listing) == 1)
        % it is not a folder ...
        fex = ~dir_listing.isdir;
    end
end
