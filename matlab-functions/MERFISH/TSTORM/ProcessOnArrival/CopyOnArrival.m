% Copy on arrival
clc;

writePath = '\\STORM4-PC\AlistairSTORM4temp\140305_tests\';
storagePath = 'T:\2014-03-05_Fluidics\';
oldFileNames = {'.';'..';'Thumbs.db'}; 

maxWait = 5; % Wait time in minutes
% Time to wait without seeing any new files appear for the system to decide its done.

t=0;
while t < maxWait
    filedata = dir(writePath);
    fileNames = {filedata.name}';

    dispatchedIdx = ismember(fileNames,oldFileNames);
    newFileNames = fileNames(~dispatchedIdx);
   
    for i=1:length(newFileNames)
    cmd = ['copy ',writePath,newFileNames{i},' ',storagePath];
       disp(cmd); 
       system(cmd); 
    end
    oldFileNames =cat(1,oldFileNames, newFileNames); 
    
    pause(60);
    t =t+1; 
    disp('watching for new files to copy...');
end