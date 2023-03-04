function RunDriftCorrect(refDaxFolder,dataDaxFolders,varargin)
% 
% 
% the 2 input arguments is a bit redundant, since the first element of
% dataDaxFolders will often be used as the reference and it is easy to
% enforce this. However 2 inputs is maybe more readable -- entry 2 will be
% registered to entry 1.  
% 
% refDaxFolder = hybFolders{1};
% dataDaxFolders = hybFolders;

defaults = cell(0,3);
defaults(end+1,:) = {'daxRoot','string','ConvZscan'}; 
defaults(end+1,:) = {'csvName','string','alignTable'};
defaults(end+1,:) = {'saveFolder','string',''};
defaults(end+1,:) = {'fidChnRef','integer',2};
defaults(end+1,:) = {'fidChnCur','integer',[]};
defaults(end+1,:) = {'overwrite','boolean',false};
defaults(end+1,:) = {'fov','integer',inf};
defaults(end+1,:) = {'LoadDaxOverwrite','boolean',false};
defaults(end+1,:) = {'LoadDaxSkipFirst','boolean',false};
pars = ParseVariableArguments(varargin,defaults,mfilename); 

if isempty(pars.fidChnCur)
    pars.fidChnCur = pars.fidChnRef;
end

if isempty(pars.saveFolder)
    analysisFolder = [refDaxFolder,'..'];
else
    analysisFolder = pars.saveFolder;
end

% find reference files
refDaxFiles = FindFiles([refDaxFolder,filesep,pars.daxRoot,'*.dax']);
nFOV = length(refDaxFiles);
nH = length(dataDaxFolders);
if isinf(pars.fov) || pars.fov==0
    fovs = 1:nFOV;
else
    fovs = pars.fov;
end

for f=fovs
    alignTableFile = [analysisFolder,pars.csvName,'_fov',num2str(f,'%03d'),'.csv'];
    if exist(alignTableFile,'file') && ~pars.overwrite
        alignTable = readtable(alignTableFile);
        if height(alignTable) >= nH
            disp(['found existing drift correction for fov',num2str(f),' using this']);
            continue
            % skip FOV that have complete alignTables already
        end
    else
        system(['del ',alignTableFile]);
    end




    [refHybFid, imProps] = LoadDax(refDaxFiles{f},'maxProject',true,'verbose',false,'channel',pars.fidChnRef,'overwrite',pars.LoadDaxOverwrite,'skipFirst',pars.LoadDaxSkipFirst); 
    for h=1:nH
        % if already exists, just load that
       runDriftCorrect = true;    
       if  exist(alignTableFile,'file')
           if height(alignTable) >= h
           fidAlign = table2struct(alignTable(alignTable.hyb==h,:));
               if ~isempty(fidAlign)
                   runDriftCorrect = false;
               end
           end
       end
               
       if runDriftCorrect
            curDaxFiles = FindFiles([dataDaxFolders{h},filesep,pars.daxRoot,'*.dax']);
            [curHybFid, imProps] = LoadDax(curDaxFiles{f},'maxProject',true,'verbose',false,'channel',pars.fidChnCur,'overwrite',pars.LoadDaxOverwrite,'skipFirst',pars.LoadDaxSkipFirst); 
            f1 = figure(1); clf; 
            fidAlign = CorrAlignFast(refHybFid,curHybFid);  
            % save data  
            fidAlign.fov = f;
            fidAlign.hyb = h;
            alignTable = struct2table(fidAlign);
            alignTableFile = [analysisFolder,pars.csvName,'_fov',num2str(f,'%03d'),'.csv'];
            writetable(alignTable,alignTableFile,'WriteMode','Append'); 
            % save image too. 
            SetFigureSavePath([analysisFolder,'CorrAlign',filesep],'makeDir',true);
            imName = ['fov',num2str(f,'%03d'),'_H',num2str(h,'%03d')];
            SaveFigure(f1,'name',imName,'formats',{'png'},'overwrite',true);
       end
    end
end     
