% AnalyzeOnArrival

dataPath = 'T:\2014-03-05_Fluidics\';
parsfile = 'T:\030414\647STORM_pars.xml';
oldDaxNames = {}; 

maxWait = 50; % Wait time in minutes
batchSize = 8;
% Time to wait without seeing any new files appear for the system to decide its done.

t=0;
while t < maxWait
    daxdata = dir([dataPath,'647*c1.dax']);
    daxNames = {daxdata.name}';
    daxFullNames = strcat(dataPath,daxNames);

    dispatchedIdx = ismember(daxNames,oldDaxNames);
    newDaxNames = daxNames(~dispatchedIdx);
    
    
    if ~isempty(newDaxNames)
        RunDotFinder('path',dataPath,'daxnames',newDaxNames,...
            'batchsize',batchSize,'parsfile',parsfile,...
            'method','DaoSTORM');
        oldDaxNames = cat(1,oldDaxNames,newDaxNames);
    end
    
    disp('watching for new files to copy...');
    pause(60);
    t =t+1; 
end