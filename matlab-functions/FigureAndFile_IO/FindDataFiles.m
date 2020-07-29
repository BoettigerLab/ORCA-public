function allFiles = FindDataFiles()
%% get list of all data files 
% this will be used to find the current drive name with data
% this is necessary because drive letters get shuffled.  

verbose = true;

drivesToSkip = {'C:','D:','E:','G:','U:','Y:'};
[~,drives] = system('wmic logicaldisk get name');
drives = strsplit(drives);
drives = drives(2:end-1);

skip = StringFind(drives,drivesToSkip);
drives(skip) = [];

allFiles = {};
for i=1:length(drives)
    filedata = dir(drives{i});
    folders = strcat(drives{i},filesep,{filedata.name})';
    if verbose
        disp(folders);
    end
    folders = folders([filedata.isdir]);
    allFiles = cat(1,allFiles{:},folders);
end




