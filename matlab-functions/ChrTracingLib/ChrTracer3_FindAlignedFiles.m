function [foundFiles,fovToSkip] = ChrTracer3_FindAlignedFiles(saveFolder,selectFOVs,numHybes,numDataChns,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true}; 
defaults(end+1,:) = {'overwrite', 'boolean', false}; 
pars = ParseVariableArguments(varargin,defaults,mfilename); 

numFOVs = length(selectFOVs);
foundFiles(numFOVs).fidNames = '';
foundFiles(numFOVs).datNames = '';
fovToSkip = false(1,numFOVs);
% ------------------- Test if data is already written ------------------
for f=selectFOVs
    fidNames = cellstr(ls([saveFolder,'fov',num2str(f,'%03d'),'*_fid.dax']));
    datNames = cell(numHybes,numDataChns);
    for n=1:numDataChns
        foundDat = cellstr(ls([saveFolder,'fov',num2str(f,'%03d'),'*_data',num2str(n),'.dax']));
        datNames(1:length(foundDat),n) = foundDat;
    end
    foundAllFid = length(fidNames) >= numHybes;
    foundAllDat = sum(cellfun(@isempty,datNames(:))) == 0;
    if foundAllFid && foundAllDat && ~pars.overwrite
        fovToSkip(f) = true;
        cprintf([1 .25 0],['found existing aligned dax files for fov ',num2str(f),' skipping alignment']);
    elseif  foundAllFid && foundAllDat && pars.overwrite
        fovToSkip(f) = false;
        cprintf([1 .25 0],['overwriting existing aligned dax files for fov ',num2str(f),' skipping alignment']);
    else
        fovToSkip(f) = false;
    end
    if ~isempty(fidNames{1})
        fidNames = strcat(saveFolder,fidNames);
    else 
        fidNames = '';
    end
    if ~isempty(datNames{1,n})
        for n=1:numDataChns
            datNames(:,n) = strcat(saveFolder,datNames(:,n));
        end
    else
       datNames = ''; 
    end
    % export FoundFiles
    foundFiles(f).fidNames = fidNames;
    foundFiles(f).datNames = datNames;
end