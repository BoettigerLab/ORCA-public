function Pipeline3(varargin)

% Global parameters 
global TSTORMdata 

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true}; 
defaults(end+1,:) = {'imageTag', 'string', 'STORM'}; 

% Parameters for DaoSTORM
defaults(end+1,:) = {'maxCPU', 'freeType', 95}; 
defaults(end+1,:) = {'batchsize', 'integer', 12}; 
defaults(end+1,:) = {'overwrite', 'integer', 0}; 
defaults(end+1,:) = {'parsfile', 'string', '\\cajal\TSTORMdata2\AdditionalAnalysis\141021_ParametersFiles\default_convpars.xml'}; 

% Thresholds to scan with DaoSTORM
defaults(end+1,:) = {'scanThresholds','boolean',true};
defaults(end+1,:) = {'thresholds','array',round(logspace(log10(200),log10(3000),20))};
defaults(end+1,:) = {'firstThreshold','positive',1000};

% Parameters for parsing file names
defaults(end+1,:) = {'fileExt', 'string', 'dax'}; % Delimiters for bin files
defaults(end+1,:) = {'fieldNames', 'cell', {'movieType', 'hybNum', 'cellNum', 'isFiducial', 'binType'}}; 
defaults(end+1,:) = {'fieldConv', 'cell', {@char, @str2num, @str2num, @(x)strcmp(x, 'c2'), @char}};
defaults(end+1,:) = {'appendExtraFields', 'bool', true};

% Figure plotting
defaults(end+1,:) = {'showFigs','string','on'};


% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

%% Main Function

thresholds = parameters.thresholds; % Threshold values to try
dataPath = parameters.dataPath; 
newPath = parameters.savePath; 


[~,initThresh] = min(abs(thresholds-parameters.firstThreshold));

%% Find all daxfiles

foundFiles = BuildFileStructure(dataPath, 'parameters', parameters,...
    'requireFlag',parameters.imageTag); % 
    % 'excludeFlags',{'WGA','HOECHST'},'requireExactMatch',true);

%% Combine dax files  
STORMmovies = strcmp({foundFiles.movieType},parameters.imageTag);
numHybes = length(unique([ foundFiles.hybNum]));
numCells = length(unique([ foundFiles.cellNum]));


% handling 
for i=1:length(foundFiles)
    if isempty(foundFiles(i).isFiducial)
        foundFiles(i).isFiducial = 0;
    end
end

targetPathInfo = dir(newPath);
if exist([dataPath,'JoinedDax'])  % && length(targetPathInfo) < 16
    robocopyCmd = ['robocopy ',[dataPath,'JoinedDax'],' ',newPath];
    system(robocopyCmd);
end

daxNames = cell(numHybes,1); 
hybeNames = cell(numHybes,1); 
for h=1:numHybes  
    daxNames{h} = ['AllCells_Hybe_',num2str(h),'.dax'];
    fileExists = exist([newPath,daxNames{h}],'file') == 2;  
    hybeN = [foundFiles.hybNum]==h-1;
    notFid = [foundFiles.isFiducial] ~= true;
    hybeNames{h} = {foundFiles.name}';
    hybeNames{h} = hybeNames{h}(hybeN & STORMmovies & notFid );
    
    if ~fileExists || parameters.overwrite == 1
        numCells = length(hybeNames{h});
        movieN = zeros(256,256,numCells); 
        for c=1:numCells
            [movieN(:,:,c),infoN] = ReadDax([dataPath,hybeNames{h}{c}],'verbose',false); 
        end
        infoN.number_of_frames = numCells;
        infoN.localPath = newPath; 
        infoN.localName = regexprep(daxNames{h},'.dax','.inf'); 
        WriteDAXFiles(movieN,infoN); 
   end
end

%% Write new parameter files with a range of threshold values

numThresholds = length(thresholds);
parsNames = cell(numThresholds,1);
for t=1:numThresholds
    fitPars = ReadParsFile(parameters.parsfile);
    fitPars.threshold = num2str(thresholds(t));
    parsNames{t} = [newPath,'convpars_t',num2str(thresholds(t)),'.xml'];
    WriteParsFile(parsNames{t} ,fitPars);
end

%% Fit dax files with the new threshold values
for t = 1:numThresholds
    % don't poll for CPU
    RunDotFinder('daxnames', strcat(newPath, daxNames),...
                 'parsfile', parsNames{t},...
                 'binname',['DAX_',num2str(thresholds(t))],...
                 'maxCPU',parameters.maxCPU,...
                 'batchsize',parameters.batchsize,...
                 'hideterminal',true,...
                 'overwrite',parameters.overwrite);
end




%% Compute the warps
disp('loading fiducial data and computing warps...'); 
% warps = cell(numCells,1);
cellProperties(numCells).warps = []; 
for c=1:numCells
    fiducialData = struct;
    fiducialData(numHybes).mList = []; 
    for h=1:numHybes
        fedNames = {foundFiles.name}'; 
        hybeN = [foundFiles.hybNum]==h-1;  % belongs to this hybe h
        isFid = [foundFiles.isFiducial] == true; % is a fiducial
        currentCell = [foundFiles.cellNum] == c-1; % belongs to current cell
        hybeName = fedNames{hybeN & STORMmovies & isFid & currentCell};
        beadBinName = regexprep(hybeName,'\.dax','_list\.bin');
        try
            fiducialData(h).mList = ReadMasterMoleculeList([dataPath,beadBinName],'verbose',false);
        catch er
            warning(er.getReport);
        end
    end
    cellProperties(c).warps = AlignFiducials(fiducialData,'verbose',false);
    % warps{c} = AlignFiducials(fiducialData,'verbose',false);
end
disp('Finished loading fiducial data and computing warps.');     


%% Load Molecule Lists, Apply Warps, Compute 
display('loading molecule lists...'); 
cellProperties(numCells).imLists = cell(numHybes,numThresholds); 
for t =1:numThresholds;
    for h = 1:numHybes
        binName = regexprep(  daxNames{h}, '\.dax',['_',num2str(thresholds(t)),'_mlist\.bin']);
        try
            imLists = ReadMasterMoleculeList([newPath,binName],'verbose',false);
        catch er
            warning(er.getReport);
        end
            
        for c=1:numCells
            cellProperties(c).imLists{h,t} = IndexStructure(imLists,imLists.frame == c);
        end
    end
end
display('Finished loading molecule lists.'); 

% copy the mfile for documentation. 
copyfile( [mfilename('fullpath'),'.m'],[newPath,mfilename,'.m']);

%% Threshold scanning!
parameters.scanThresholds = true;

try
    [bestThetas,initRho,finalRho] = ScanThresholds2(cellProperties,thresholds,initThresh,newPath,'parameters',parameters);
    disp(['initial FPKM correlation = ',num2str(initRho)]);
    disp(['final FPKM correlation = ',num2str(finalRho)]);
catch er
    warning(er.getReport); 
    parameters.scanThresholds = true;
    [bestThetas,initRho,finalRho] = ScanThresholds2(cellProperties,thresholds,initThresh,newPath,'parameters',parameters);
    disp(['initial FPKM correlation = ',num2str(initRho)]);
    disp(['final FPKM correlation = ',num2str(finalRho)]);
end


%% Save DATA
% -------------------------------------------------------------

% save molecule lists
numCells = length(cellProperties);
thetas = bestThetas;%  5*ones(1,numHybes);
for c=1:numCells
    idx = sub2ind([numHybes,numThresholds], (1:numHybes)',thetas'); 
    imLists = cellProperties(c).imLists(idx);
    for h=1:numHybes
        newMListName = regexprep(hybeNames{h}{c},'\.dax','_final_alist\.bin');
        imLists{h}.frame = ones(length(imLists{h}.xc),1,'int32'); % reset all frames to frame 1
        WriteMoleculeList(imLists{h},[dataPath,newMListName],'verbose',false); 
        WriteMoleculeList(imLists{h},[newPath,newMListName],'verbose',false); 
    end
      disp(['finished writing molecule lists for cell ',num2str(c),' of ',num2str(numCells)]);
end    
