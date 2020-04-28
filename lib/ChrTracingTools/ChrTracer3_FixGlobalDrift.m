function [fiducialAlignFrames,regData] = ChrTracer3_FixGlobalDrift(fiducialFrames,varargin)


% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
global scratchPath;
defaults = cell(0,3);
% key general parameters
defaults(end+1,:) = {'selectFOVs','freeType',inf}; 
defaults(end+1,:) = {'saveFits', 'boolean', true}; 
defaults(end+1,:) = {'saveFolder', 'string', scratchPath}; 
% key alignment parameters
defaults(end+1,:) = {'alignContrastLow', 'fraction', .7}; % low image threshold for contrast balance prior to coarse alignment
defaults(end+1,:) = {'alignContrastHigh', 'fraction', .9995}; % high threshold  for contrast balance prior to coarse alignment
defaults(end+1,:) = {'refHybe','integer',1}; % hybe to use to start alignment
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

%=== main function
[numHybes,numFOVs] = size(fiducialFrames);
fiducialAlignFrames = cell(numFOVs,1);
regData = cell(numFOVs,1);


selectFOVs = true(1,numFOVs);
if any(pars.selectFOVs) && ~any(isinf(pars.selectFOVs))
    selectFOVs  = ~selectFOVs;
    selectFOVs(pars.selectFOVs) = true;
end

%---------- Check for existing data, skip if necessary
%=== becomes obsolete in ChrTracer3p2
% skip FOV that already have aligned data files
[~,fovToSkip] = ChrTracer3_FindAlignedFiles(pars.saveFolder,selectFOVs,numHybes,1);
fovToSkip(pars.skipFOVs) = true;
if ~pars.overwrite
    if pars.verbose && any(fovToSkip)
       disp(['Found existing aligned files, skipping fov ',num2str(find(fovToSkip))]); 
    end
    selectFOVs(fovToSkip) = false; % remove previously analyzed fovs from list of to drift correct   
end

% ask if we should skip FOV that already have _regData.csv files
alsoSkip = false(1,length(selectFOVs));
for f=find(selectFOVs)
    tableName = [pars.saveFolder,'fov',num2str(f,'%03d'),'_regData.csv'];
    if exist(tableName,'file') && ~pars.overwrite
       answer = questdlg(['found a saved alignment map for FOV ',num2str(f),'.  load it and skip this FOV?'], ... % question
                        'Skip', ...  % pop-up label
                        'Yes','No','Yes'); % op1 op2 default
       if strcmp(answer,'Yes')
           alsoSkip(f) = true;
           tableData = readtable(tableName);
           regData{f} = table2struct(tableData);
       end
    end
end
selectFOVs(alsoSkip) = false;

% --- actually registered the frames and save a .csv of the result
%   (in future we should try just using the csv file and not rewriting the
%   data as an aligned copy). 
for f=find(selectFOVs)
    [fiducialAlignFrames{f},regData{f}] = RegisterImagesFast(fiducialFrames(:,f),'parameters',pars);
    % save data
    if pars.saveFits
        tableData = struct2table(regData{f});
        tableName = [pars.saveFolder,'fov',num2str(f,'%03d'),'_regData.csv'];
        writetable(tableData,tableName);
    end
end

