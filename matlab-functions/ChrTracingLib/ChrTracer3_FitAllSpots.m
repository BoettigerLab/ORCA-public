function spotDataTable = ChrTracer3_FitAllSpots(fidSpots,datSpots,varargin)

defaults = cell(0,3);
% FitAllSpots default parameters
defaults(end+1,:) = {'numParallel','integer',1};
defaults(end+1,:) = {'selectSpots','integer',[]};
defaults(end+1,:) = {'saveTable','boolean',true};

% Parameters forwarded to ChrTracer3_FitSpots
% General FOV parameters
defaults(end+1,:) = {'nmXYpix','positive',154};
defaults(end+1,:) = {'nmZpix','positive',100};
defaults(end+1,:) = {'fov','nonnegative',1};
defaults(end+1,:) = {'lociXY', 'array', [0,0]};
defaults(end+1,:) = {'saveFolder','string',''};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'showPlots','boolean',false};
defaults(end+1,:) = {'overwrite','boolean',false};
defaults(end+1,:) = {'saveData','boolean',false}; % changed
defaults(end+1,:) = {'eTable','freeType',[]};

% Fiducial alignment defaults
defaults(end+1,:) = {'upsample','positive',4};  % 8 for accuracy 2 for speed
defaults(end+1,:) = {'maxXYdrift','positive',4};
defaults(end+1,:) = {'maxZdrift','positive',6};
defaults(end+1,:) = {'fidRegTheta','nonnegative',.6};

% Fiducial fitting defaults
defaults(end+1,:) = {'fidMinPeakHeight', 'positive', 200};
defaults(end+1,:) = {'fidCameraBackground', 'nonnegative', 0};
defaults(end+1,:) = {'fidPeakBlur', 'nonnegative', .5};
defaults(end+1,:) = {'fidTroubleshoot', 'boolean', false};
defaults(end+1,:) = {'fidMaxFitWidth', 'positive', 8}; % 6
defaults(end+1,:) = {'fidMaxFitZdepth', 'positive', 12}; % 6
defaults(end+1,:) = {'fidMinSep', 'nonnegative', 5};  % Min separation between peaks in pixels.  Closer than this will be averaged
defaults(end+1,:) = {'fidKeepBrightest','integer',2};  % Max number of peaks to allow
defaults(end+1,:) = {'fidRelativeHeight','fraction',0};

% Data fitting defaults
defaults(end+1,:) = {'nmXYpix','positive',154};
defaults(end+1,:) = {'nmZpix','positive',100};
defaults(end+1,:) = {'maxXYstep','positive',5}; % box radius in pixels
defaults(end+1,:) = {'maxZstep','positive',7}; % box radius in pixels
defaults(end+1,:) = {'datMinPeakHeight', 'positive', 500};
defaults(end+1,:) = {'datCameraBackground', 'nonnegative', 0};
defaults(end+1,:) = {'datPeakBlur', 'nonnegative', .5};
defaults(end+1,:) = {'datTroubleshoot', 'boolean', false};
defaults(end+1,:) = {'datMaxFitWidth', 'positive', 8};
defaults(end+1,:) = {'datMaxFitZdepth', 'positive', 12};
defaults(end+1,:) = {'datMinSep', 'nonnegative', 1};
defaults(end+1,:) = {'datMinHBratio','nonnegative',1.2}; % peak value over background value
defaults(end+1,:) = {'datMinAHratio','nonnegative',.25}; % fitted height over background vs peak value
defaults(end+1,:) = {'datMaxUncert','nonnegative',2}; % pixels
defaults(end+1,:) = {'keepBrightest','integer',1}; % will use variable idx to return the brightness rank of each spot  

% for debuggin
defaults(end+1,:) = {'stopOnError','boolean',true};


% Chromatic Align parameters
defaults(end+1,:) = {'refChn','integer',647};
defaults(end+1,:) = {'minHeightRef','nonnegative',0};
defaults(end+1,:) = {'minHeightOt','nonnegative',0};
% parameters for Polymap3D
defaults(end+1,:) = {'chnColors','colormap',hsv(3)};
defaults(end+1,:) = {'polyOrder','integer',2};
% defaults(end+1,:) = {'showPlots','boolean',true};  % this par is global

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



tic
dataTable = cell(numSpots,1);
% select a spot and extract the image data 
f = pars.fov;
if pars.numParallel < 2
    for s = selectSpots % 320;
        fidSpts = fidSpots{s};
        datSpts = datSpots{s};
        try
            dataTable{s} = ChrTracer3_FitSpots(fidSpts,datSpts,s,'parameters',pars);
        catch er
            if pars.verbose
                warning(['encountered error on fov ',num2str(f), ' spot ',num2str(s)]);
                warning(er.getReport);
            end
            if pars.stopOnError
                disp('set breakpoint here');              
            end
        end
    end
else
    parfor s=selectSpots
        fidSpts = fidSpots{s};
        datSpts = datSpots{s};
        try
        dataTable{s} = ChrTracer3_FitSpots(fidSpts,datSpts,s,'parameters',pars);
        catch 
           if pars.verbose
                warning(['encountered error on fov ',num2str(f), ' spot ',num2str(s)]);
           end 
        end
    end
end

% combine all the tables
spotDataTable = cat(1,dataTable{:});

% compute chromatic warp from Align Data for this FOV. 
%   Note, it is better to combine all data from multiple FOVs and use that.
nChns = size(datSpots{1},2);
if nChns > 1
    try
    spotDataTable = ChromaticAlignFromTable(spotDataTable,'parameters',pars);
    catch er
        disp(er.message);
        disp(['WARNING: ERROR running ChromaticAlignFromTable for fov ',num2str(pars.fov,'%03d')]);
    end
end

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
