function filespath = getdatesortedfiles(directory)
%This function returns the latest file from the directory passsed as input
%argument

%Get the directory contents
dirc = dir(directory);

fileListArray = cellfun(@(c)[directory c], {dirc(:).name}, 'uni',false);

%Filter out all the folders.
dirc = dirc(find(~cellfun(@isdir,fileListArray)));

%I contains the index to the biggest number which is the latest file
[A,I] = sort([dirc(:).datenum]);

filespath = {};
if ~isempty(I)    
    for (i = 1:length(I))
        filespath{end+1} = [dirc(i).folder, '/', dirc(i).name];
    end
end

end