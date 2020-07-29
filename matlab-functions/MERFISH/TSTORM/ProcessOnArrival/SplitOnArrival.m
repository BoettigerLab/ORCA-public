
storagePath = 'L:\140310_L4E4\latecopy\';
oldDaxNames = {}; 

maxWait = 15; % Wait time in minutes
% Time to wait without seeing any new files appear for the system to decide its done.

t=0;
while t < maxWait
    daxdata = dir([storagePath,'*.dax']);
    daxNames = {daxdata.name}';
    daxFullNames = strcat(storagePath,daxNames);

    dispatchedIdx = ismember(daxNames,oldDaxNames);
    newDaxNames = daxdata(~dispatchedIdx);
    
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
    
    disp('watching for new files to split...');
    pause(60);
    t =t+1; 
end