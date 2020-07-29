% Copy and analyze on arrival
clc;

writePath = '\\STORM4-PC\AlistairSTORM4temp\140317_L3E7\'; % where the raw data appears / is first saved
storagePath = 'L:\140317_L3E7\'; % where the raw data should be copied too
oldFileNames = {'.';'..';'Thumbs.db'}; 

dataPath = 'L:\140317_L3E7\';
parsfile = 'L:\140317_L3E7\647iniPars.ini';
parsfile2 ='L:\140317_L3E7\AveBeadPars.ini';
oldDaxNames = {}; 
oldDaxNames2 = {}; 

maxWait = 60; % Wait time in minutes
batchSize = 4;

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
    
    % now look of for dax files with c1 label
    daxdata = dir([dataPath,'*c1.dax']);
    daxNames = {daxdata.name}';

    dispatchedIdx = ismember(daxNames,oldDaxNames);
    newDaxNames = daxNames(~dispatchedIdx);
    
    if ~isempty(newDaxNames)
        RunDotFinder('path',dataPath,'daxnames',newDaxNames,...
            'batchsize',batchSize,'parsfile',parsfile,...
            'method','insight','overwrite',0,'hideterminal',true);
        oldDaxNames = cat(1,oldDaxNames,newDaxNames);
    end
    
    % now look of for dax files with c2 label
    daxdata2 = dir([dataPath,'*c2.dax']);
    daxNames2 = {daxdata2.name}';

    dispatchedIdx = ismember(daxNames2,oldDaxNames2);
    newDaxNames2 = daxNames2(~dispatchedIdx);
    
    if ~isempty(newDaxNames2)
        RunDotFinder('path',dataPath,'daxnames',newDaxNames2,...
            'batchsize',batchSize,'parsfile',parsfile2,...
            'method','insight','overwrite',0,'hideterminal',true);
        oldDaxNames2 = cat(1,oldDaxNames2,newDaxNames2);
    end
    
    disp(['watching for new files to copy or analyze...',' t=',num2str(t)]);
    pause(60);
    t =t+1; 
end