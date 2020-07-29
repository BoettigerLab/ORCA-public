function [fiducialAlignFrames,regData] = ChrTracer3p2_FixGlobalDrift(fiducialFrames,varargin)


% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
global scratchPath;
defaults = cell(0,3);
% key general parameters
defaults(end+1,:) = {'selectFOVs','freeType',inf}; 
defaults(end+1,:) = {'saveFits', 'boolean', true}; 
defaults(end+1,:) = {'saveFigs', 'boolean', true}; 
defaults(end+1,:) = {'loadPriorData', 'boolean', true}; 
defaults(end+1,:) = {'saveFolder', 'string', scratchPath}; 
% key alignment parameters
defaults(end+1,:) = {'alignContrastLow', 'fraction', .7}; % low image threshold for contrast balance prior to coarse alignment
defaults(end+1,:) = {'alignContrastHigh', 'fraction', .9995}; % high threshold  for contrast balance prior to coarse alignment
defaults(end+1,:) = {'refHybe','integer',1}; % hybe to use to start alignment
defaults(end+1,:) = {'hybs','integer',inf}; % hybe to use to start alignment
% other CorrAlignFast Parameters
defaults(end+1,:) = {'maxSize', 'positive', 400}; % rescale all images to this size for alignment
defaults(end+1,:) = {'fineBox', 'freeType', 100};  % perform fine scale alignment using a box of this size around the brightest point.
defaults(end+1,:) = {'fineUpsample', 'positive', 3};  
defaults(end+1,:) = {'maxShift', 'nonnegative', inf};
defaults(end+1,:) = {'gradMax', 'boolean', true};
defaults(end+1,:) = {'minGrad', 'float', -inf};
defaults(end+1,:) = {'angles','float',0}; % -10:1:10
defaults(end+1,:) = {'scales','float',1}; % -10:1:10
defaults(end+1,:) = {'fineMaxShift', 'nonnegative', 10};
defaults(end+1,:) = {'fineAngles','float',0}; % -1:.1:1
defaults(end+1,:) = {'fineScales','float',1}; % 0.95:0.01:.1.05
defaults(end+1,:) = {'fineCenter','array',[0,0]};
defaults(end+1,:) = {'showplot', 'boolean', true};
defaults(end+1,:) = {'fastDisplay', 'boolean', true};
defaults(end+1,:) = {'displayWidth', 'integer', 500};
defaults(end+1,:) = {'showExtraPlot', 'boolean', false};
defaults(end+1,:) = {'minFineImprovement', 'float', 0}; 
defaults(end+1,:) = {'showCorrAlign', 'boolean', true};
% common fov parameters
defaults(end+1,:) = {'overwrite', 'boolean', false}; 
defaults(end+1,:) = {'verbose', 'boolean', true}; 
defaults(end+1,:) = {'veryVerbose', 'boolean', false}; 
defaults(end+1,:) = {'showPlots', 'boolean', true}; 
defaults(end+1,:) = {'showExtraPlots', 'boolean', false}; 
defaults(end+1,:) = {'stopOnError','boolean',false};
% internal
defaults(end+1,:) = {'skipFOVs','freeType',[]}; 
pars = ParseVariableArguments(varargin,defaults,mfilename);

% enforce saveFolder must end with a filesep.
if ~isempty(pars.saveFolder)
   if ~strcmp(pars.saveFolder(end),filesep)
      pars.saveFolder = [pars.saveFolder,filesep]; 
   end
end

if pars.saveFigs
    SetFigureSavePath(pars.saveFolder,'verbose',false);
end

% note, size of fiduciaFrames is determined by the number of .dax files
% matching the daxroot (default ConvZscan) in the file folder. 
% if a subset of FOVs is chosen using the FOV select, these should be stuck
% into the corresponding positions of the fiducialFrames, leaving blank
% cells anywhere data was skipped. 

%=== main function
[numHybes,numFOVs] = size(fiducialFrames); % 
fiducialAlignFrames = cell(numFOVs,1);
regData = cell(numFOVs,1);

% search for existing files
foundFiles = cellstr(ls([pars.saveFolder,'*_regData.csv']));
hasData = false(1,numFOVs);
if ~isempty(foundFiles{1})   
    idHasData = cellfun(@(x) str2num(x(4:6)),foundFiles);
    hasData(idHasData) = true; 
end

% handle select FOVs
selectFOVs = true(1,numFOVs);
if any(pars.selectFOVs) && ~any(isinf(pars.selectFOVs))
    selectFOVs  = ~selectFOVs;
    selectFOVs(pars.selectFOVs) = true;
end

% ask if redo?
answer = 'Overwrite All'; % if there is nothing to overwrite, there is no prompt 
needRedo = hasData & selectFOVs;
idxRedo = find(needRedo)
if sum(needRedo) > 0
   disp('found existing registration files for FOVs: ');
   disp(idxRedo);
   disp('select how to handle existing files. Your selection will be applied after all overlaps are resolved');
   answer = questdlg(['Registration files exist'], ... % question
                            'Skip', ...  % pop-up label
                            'Overwrite All','Skip All','Overwrite Selected','Skip All'); % op1 op2 default.  max of 3 buttons!
end

alsoSkip = false(1,length(selectFOVs));
if strcmp(answer,'Skip All')
    for f=idxRedo   % need to load data for any FOV with data
        tableName = [pars.saveFolder,'fov',num2str(f,'%03d'),'_regData.csv'];
        tableData = readtable(tableName);
        regData{f} = table2struct(tableData); % record data
        alsoSkip(f) = true;  % no need to align
    end
    selectFOVs(alsoSkip) = false;
elseif strcmp(answer,'Overwrite All')
    disp('will now overwrite data registration data for FOVs');
    disp(idxRedo);
elseif strcmp(answer,'Overwrite Selected')
    toOverwrite = input('which FOVs #s should be overwritten (enter in "quotes")? ');
    toOverwrite = str2num(toOverwrite); %#ok<ST2NM>
    toLoad = needRedo;
    toLoad(toOverwrite) = false;
    for f=find(toLoad) % load anything that was not going to be overwritten
        tableName = [pars.saveFolder,'fov',num2str(f,'%03d'),'_regData.csv'];
        tableData = readtable(tableName);
        regData{f} = table2struct(tableData); % record data
        alsoSkip(f) = true;  % no need to align
    end
    selectFOVs(alsoSkip) = false; % anything we didn't load will rerun
end


% --- actually registered the frames and save a .csv of the result
%   (in future we should try just using the csv file and not rewriting the
%   data as an aligned copy). 
for f=find(selectFOVs)
    disp(['aligning FOV ',num2str(f)]);
    if isempty(fiducialFrames(1,f))
        error('fiducialFrames is empty, did you skip loading max projection data in step 1?');
    end
    [fiducialAlignFrames{f},regData{f}] = RegisterImagesFast(fiducialFrames(:,f),'parameters',pars,'fov',f);
    % save data
    if pars.saveFits
        tableData = struct2table(regData{f}); % annyoingly, this groups  
        tableName = [pars.saveFolder,'fov',num2str(f,'%03d'),'_regData.csv'];
        if exist(tableName,'file')==2 % just update those Hybs which ran alignment
            currTable = readtable(tableName);
            notUpdated = isnan(tableData{:,1});
            tableDataOut = tableData;
            tableDataOut{notUpdated,:} = currTable{notUpdated,:};
        else
            tableDataOut = tableData;
        end
        
        writetable(tableDataOut,tableName); % will overwrite
    end
    if pars.saveFigs
       SaveFidAlignFrames(fiducialAlignFrames{f},'fov',f,'overwrite',true,'verbose',pars.verbose);
       % if overwriting the table, we should overwrite the image
    end
end

if pars.verbose
    disp(['finished aligning ', num2str(sum(selectFOVs)), ' FOVs']);
end

