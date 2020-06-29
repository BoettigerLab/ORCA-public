function spotDataTable = ChrTracer2_FitAllSpots3(fidSpots,datSpots,varargin)

defaults = cell(0,3);
% FitAllSpots default parameters
defaults(end+1,:) = {'numParallel','integer',1};
defaults(end+1,:) = {'selectSpots','integer',[]};
defaults(end+1,:) = {'saveTable','boolean',true};

% Parameters forwarded to ChrTracer2_FitSpots
% General FOV parameters
defaults(end+1,:) = {'nmXYpix','positive',154};
defaults(end+1,:) = {'nmZpix','positive',100};
defaults(end+1,:) = {'fov','nonnegative',1};
defaults(end+1,:) = {'lociXY', 'array', [0,0]};
defaults(end+1,:) = {'saveFolder','string',''};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'showPlots','boolean',false};
defaults(end+1,:) = {'overwrite','boolean',false};
defaults(end+1,:) = {'eTable','freeType',[]};
defaults(end+1,:) = {'saveData','boolean',true};

% Fiducial alignment defaults
defaults(end+1,:) = {'upsample','positive',4};  % 8 for accuracy 2 for speed
defaults(end+1,:) = {'maxShiftXY','positive',4};
defaults(end+1,:) = {'maxShiftZ','positive',6};
defaults(end+1,:) = {'fidRegTheta','nonnegative',.6};

% Fiducial fitting defaults
defaults(end+1,:) = {'fidMinPeakHeight', 'positive', 200};
defaults(end+1,:) = {'fidCameraBackground', 'nonnegative', 0};
defaults(end+1,:) = {'fidPeakBlur', 'nonnegative', .5};
defaults(end+1,:) = {'fidTroubleshoot', 'boolean', false};
defaults(end+1,:) = {'fidMaxFitWidth', 'positive', 8};
defaults(end+1,:) = {'fidMinSep', 'nonnegative', 5};  % Min separation between peaks in pixels.  Closer than this will be averaged
defaults(end+1,:) = {'fidKeepBrightest','integer',1};  % Max number of peaks to allow
defaults(end+1,:) = {'fidRelativeHeight','fraction',0};

% Data fitting defaults
defaults(end+1,:) = {'datBoxXY','positive',5}; % box radius in pixels
defaults(end+1,:) = {'datBoxZ','positive',7}; % box radius in pixels
defaults(end+1,:) = {'datMinPeakHeight', 'positive', 500};
defaults(end+1,:) = {'datCameraBackground', 'nonnegative', 0};
defaults(end+1,:) = {'datPeakBlur', 'nonnegative', .5};
defaults(end+1,:) = {'datTroubleshoot', 'boolean', false};
defaults(end+1,:) = {'datMaxFitWidth', 'positive', 6};
defaults(end+1,:) = {'datMinSep', 'nonnegative', 3};
defaults(end+1,:) = {'datMinHBratio','nonnegative',1.2}; % peak value over background value
defaults(end+1,:) = {'datMinAHratio','nonnegative',.25}; % fitted height over background vs peak value
defaults(end+1,:) = {'datMaxUncert','nonnegative',2}; % pixels


% parse defaults
pars = ParseVariableArguments(varargin,defaults,mfilename); 
% pars = ParseVariableArguments([],defaults,'FitSpots2'); 
%%

disp('Fitting Data...');

% % TEST MODE loading test data straight from CT
% global CT
% fidSpots = CT{1}.fidSpots;
% datSpots = CT{1}.dataSpots;
% pars.lociXY = CT{1}.parsFOV.lociXY;
% pars.fov = CT{1}.parsFOV.fov;
% pars.saveFolder = CT{1}.parsFOV.saveFolder;

numSpots = length(fidSpots);
if isempty(pars.selectSpots)
    selectSpots = 1:numSpots;
else
    selectSpots = pars.selectSpots;
end


if pars.numParallel > 1
   p = gcp('nocreate');
   if isempty(p)
      parpool(pars.numParallel); 
      p = gcp('nocreate');
   elseif pars.numParallel ~= p.NumWorkers
      delete(gcp('nocreate')); 
      parpool(pars.numParallel); 
      p = gcp('nocreate');
   end
   if pars.verbose
      disp(['Using ',num2str(p.NumWorkers),' cores.']); 
   end
end


tic
dataTable = cell(numSpots,1);
% select a spot and extract the image data 
if pars.numParallel < 2
    for s = selectSpots % 320;
        fidSpts = fidSpots{s};
        datSpts = datSpots{s,1};
        dataTable{s} = ChrTracer2_FitSpots3(fidSpts,datSpts,s,'parameters',pars,'veryverbose',true);
    end
else
    parfor s=selectSpots
        fidSpts = fidSpots{s};
        datSpts = datSpots{s,1};
        dataTable{s} = ChrTracer2_FitSpots3(fidSpts,datSpts,s,'parameters',pars);
    end
end

spotDataTable = cat(1,dataTable{:});
% save table 
if pars.saveData && pars.saveTable
    tableSaveName = [pars.saveFolder, 'fov',num2str(pars.fov,'%03d'),'_AllFits.csv'];
    writetable(spotDataTable,tableSaveName);
    disp(['wrote ',tableSaveName]);
end
t = toc;

if pars.verbose
    disp(['FitSpots completed in ',num2str(t/60),' minutes']);
end

%%
