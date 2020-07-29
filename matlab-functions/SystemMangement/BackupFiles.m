% Back stuff up
%% 
dataDrive = 'L:\';
backupDrive = 'F:\RAID1backup';

excludeDax =  ' /xf *.dax';
copySubfolders = ' /s';
excludeSTORM = ' /max:101000000'; % 101000000  250000000 - 24MB   1200000 1 % 
retries = ' /r:1 /w:2'; 
notHidden = ' /A-:HASRNT'; % turn off Hidden, Archive, System, Read only, Not content indexed, Temporary  
command = ['robocopy ',dataDrive,' ',backupDrive, excludeSTORM,copySubfolders,retries,notHidden,' &'];
disp(command);
system(command);
system(['attrib -h -s ',backupDrive]); % make folder visible. 
% Otherwise Robocopy makes this invisible after copying.  


%% Old approach
% 
% backupDrive = 'I:';
% searchFolders = {...
%     'T:',...
%     }
% 
% i= 1;
% % find subfolders
% allDir = dir(searchFolders{i});
% 
% subFolders = {allDir.name};
% subFolders(~[allDir.isdir]) = [];
% 
% 
% for s = 1:length(subFolders) % s=6
%     searchFolder = [searchFolders{i},filesep,subFolders{s},filesep];
%     listSubFolders = genpath(searchFolder);
%     dataFolders = strsplit(listSubFolders,';')';
%     
%     for ss=1:length(dataFolders)  % ss = 3
%         alldata = dir(dataFolders{ss});
%         isDax = ~cellfun(@isempty, strfind({alldata.name},'.dax'));
%         dataSize = [alldata.bytes];
%         sum(dataSize(~isDax))/sum(dataSize)
%         
%         
%         allNames = {alldata.name};
%         isSpecial = strcmp(allNames,'..');
%         allNames( [alldata.isdir] | isDax | isSpecial) = [];
%         
%         for f = 1:length(allNames)
%             oldFolder = dataFolders{ss};
%             newFolder = regexprep(dataFolders{ss},searchFolders{i},backupDrive);
%             origFile =  [oldFolder,filesep,allNames{f}];
%             newFile =   [newFolder,filesep,allNames{f}];
%         end
%         
%         
%     end
% end
