function LoadMiji

warning('All global variables will be erased!');
cont = input('continue? (1=yes,0=no)  ');

if cont
% javaaddpath 'C:\Program Files\MATLAB\R2014a\java\mij.jar'.
    global  mijiPath; 
    javaaddpath(mijiPath);  % this clears all paths 
    ResetPaths;
    global fijiPath; % this clears all paths
    addpath(fijiPath); 
    Miji(0);  % Miji seems to wipe my paths
    ResetPaths; 
end