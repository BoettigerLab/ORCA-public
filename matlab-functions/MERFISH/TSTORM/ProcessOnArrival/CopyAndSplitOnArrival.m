% Copy and split on arrival
clc;

writePath = '\\STORM4-PC\AlistairSTORM4temp\140313_L3E8\'; % where the raw data appears / is first saved
storagePath = 'L:\140313_L3E8\'; % where the raw data should be copied too
oldFileNames = {'.';'..';'Thumbs.db'}; 

oldDaxNames = {}; 
maxWait = 60; % Wait time in minutes


% Time to wait without seeing any new files appear for the system to decide its done.

t=0;
while t < maxWait
    filedata = dir(writePath);
    fileNames = {filedata.name}';
    fileSizes = [filedata.bytes]';
    fileNames = fileNames(fileSizes > 10);

    dispatchedIdx = ismember(fileNames,oldFileNames);
    newFileNames = fileNames(~dispatchedIdx);
   
    for i=1:length(newFileNames)
        if ~exist([storagePath,newFileNames{i}],'file')
        cmd = ['copy ',writePath,newFileNames{i},' ',storagePath];
        disp(cmd); 
        system(cmd); 
        else 
            disp(['File ',newFileNames{i},' already exists in target folder. Skipping']); 
        end
        t=0;
    end
    oldFileNames =cat(1,oldFileNames, newFileNames); 
    
    daxdata = dir([storagePath,'*.dax']);
    daxNames = {daxdata.name}';
    daxFullNames = strcat(storagePath,daxNames);

    dispatchedIdx = ismember(daxNames,oldDaxNames);
    newDaxNames = daxNames(~dispatchedIdx);
    %newDaxNames = daxdata(~dispatchedIdx);
    
    savePath = [storagePath,'splitdax\'];
    if ~exist(savePath,'dir')
        mkdir(savePath);
    end
    
    if ~isempty(newDaxNames)
        splitQVdax(storagePath,'alldax',newDaxNames,...
            'savepath',savePath,...
          'chns',{'647','561'},'delete',true,'minsize',5E2);
        oldDaxNames = cat(1,oldDaxNames,newDaxNames);
    end
    
    disp(['watching for new files to copy or analyze...',' t=',num2str(t)]);
    pause(60);
    t =t+1; 
end