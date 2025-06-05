%% Fit all RNA files in folder
% clear all; close all; startup;

expFolder = '\\169.254.113.81\NAS02_Vol1\Sedona\2023-08-16_96-168hrgastruloids_hoxRNA_try2\';
rnaFolder = [expFolder, 'RNA_Expt\'];  
rnaTable = readtable([rnaFolder,'Hyb_Identities.xlsx'])  
rnaNames = rnaTable.RNA

% hybFolders = FindFiles([rnaFolder,'Hyb*'],'onlyFolders',true); % autofind 
hybFolders = strcat(rnaFolder,rnaTable.FolderName,'/') ; % in this case Hybs 1-5 are membrane barcodes, not RNA
hybFolders = hybFolders(1:16) % only the first 16 are RNA data
nH = length(hybFolders);

analysisFolder = SetFigureSavePath([rnaFolder,'Analysis\'],'makeDir',true);
rnaFitFolder = [analysisFolder,'RNA_Fits\']; % RNA fit data will be saved here
nH = length(hybFolders);
daxFiles_h1 = FindFiles([hybFolders{1},'*.dax']);
nFOV = length(daxFiles_h1);
if nFOV ==0
    error(['no dax files found in folder ',hybFolders{1}]);
end
numParallel = 20;


%% find all orig daxFiles
daxFiles = cell(nFOV,nH);
for h=1:nH
    temp = FindFiles([hybFolders{h},'C*.dax']);
    daxFiles(1:length(temp),h) = temp;
end
%%
maxProj = false;
r=0;
runOut = cell(1,1); 
verbose = false;
runSilent = true; 
for f=1:nFOV
    k=0;
    for h=1:nH
        disp(['processing FOV ',num2str(f),' of ',num2str(nFOV), '  hyb ',num2str(h), ' of ',num2str(nH), ' ...']);
        k=k+1;
        saveFolder = [analysisFolder,'RNA_Fits\'];
        
        daxNewName1 = ['3D_Hyb',rnaNames{h,1},'_FOV',num2str(f,'%02d')]; % full 3D
      
        daxNewFile1 = [saveFolder,daxNewName1,'.dax'];
        binName = regexprep(daxNewFile1,'.dax','.hdf5');
        csvName = regexprep(daxNewFile1,'.dax','.csv');
        
        if  ~exist(csvName,'file') % ~exist(binName,'file') &&
            if exist(binName,'file')
                disp('rewriting binfile')
                delete(binName);
                pause(.1);
            end

            try
                [daxIn,daxInfo1] = ReadDax(daxFiles{f,h},'verbose',verbose);
                daxInfo2 = daxInfo1;
            catch er
                warning(['Failed to load dax FOV=',num2str(f),'  hyb=',num2str(h),'  name=',daxFiles{f,h}]);
                % disp(er.getReport);
                continue;
            end

            daxOut1 = daxIn(:,:,1:2:end); % this assumes 2-channel, alternating data
            daxOut1 = daxOut1(:,:,3:end-2); % toss out start and end frames;
            WriteDax(daxOut1, 'daxName',daxNewName1,'infoFile',daxInfo1,'folder',saveFolder,'confirmOverwrite',false,'verbose',verbose);
            pause(1); 
            r=r+1;
            if ~exist([saveFolder,'tempDaoPars.xml'],'file')
                runOut{r} = DaoFit(daxNewFile1,'overwrite','resume','max_frame',-1,...
                       'threshold',20,'camera_offset',100,'sigma',1.3,...
                       'foreground_sigma',1,'background_sigma',10,...
                       'verbose',verbose,'showPlots',false,'runExternal',true,'Hidden',runSilent);
            else
                runOut{r} = DaoFit(daxNewFile1,'parsFile',[saveFolder,'tempDaoPars.xml'],...
                       'verbose',verbose,'showPlots',false,'runExternal',true,'Hidden',runSilent);
            end
            isRunning = ~cellfun(@(x) x.HasExited,runOut);
            numRunning = sum(isRunning);
            while numRunning > numParallel
                pause(1);
                isRunning = ~cellfun(@(x) x.HasExited,runOut);
                numRunning = sum(isRunning);
            end
        end
    end


    %% convert
    binNames = FindFiles([rnaFitFolder,'*.hdf5']);
    nBins = length(binNames);
    r = 0;
    runOut = {};
    for b=1:nBins
        csvFile = regexprep(binNames{b},'.hdf5','.csv');
        % if a csv-version of the data has not already been created, generate
        % it using the storm_analysis converter:
        if ~exist(csvFile,'file')
            pythonSA = [condaPrompt,'cd ', pyPath, 'Lib\site-packages\storm_analysis\ && python.exe '];
            hd5_to_txt = './sa_utilities/hdf5_to_txt.py';
            cmdOut = [pythonSA,hd5_to_txt ' --hdf5 ',binNames{b},' --txt ',csvFile];
            r=r+1;
            runOut{r} = SystemRun(cmdOut,'Hidden',true);
            isRunning = ~cellfun(@(x) x.HasExited,runOut);
            numRunning = sum(isRunning);
            while numRunning > numParallel
                pause(1);
                isRunning = ~cellfun(@(x) x.HasExited,runOut);
                numRunning = sum(isRunning);
            end
        end
    end

    if ~isempty(runOut)
        if ~isempty(runOut{1})
            % wait for all processes to complete before running clean-up
           isRunning = ~cellfun(@(x) x.HasExited,runOut);
            numRunning = sum(isRunning);
            while numRunning >= 1
                isRunning = ~cellfun(@(x) x.HasExited,runOut);
                numRunning = sum(isRunning);
                disp(['waiting for fits to complete...'])
                pause(10);
            end   
            % ---- cleanup
            daxNames = FindFiles([rnaFitFolder,'*.dax']);
            for d=1:length(daxNames)
                cmd = ['del ',daxNames{d}];
                system(cmd);
            end   
            daxNames = FindFiles([rnaFitFolder,'*.dax']);
        end
    end
end

%% remove all temporary dax and hdf5 files

daxNames = FindFiles([rnaFitFolder,'*.dax'])
for d=1:length(daxNames)
    cmd = ['del ',daxNames{d}];
    system(cmd);
end

hd5Names = FindFiles([rnaFitFolder,'*.hdf5']);
for d=1:length(hd5Names)
    cmd = ['del ',hd5Names{d}];
    system(cmd);
end
