function spotDataTable = ChrTracer3_MultiFitAllSpots(fidSpots,datSpots,varargin)

defaults = cell(0,3);
% FitAllSpots default parameters
defaults(end+1,:) = {'numParallel','integer',1};
defaults(end+1,:) = {'selectSpots','integer',[]};
defaults(end+1,:) = {'saveTable','boolean',true};


% General FOV parameters
defaults(end+1,:) = {'fov','nonnegative',1};  % key
defaults(end+1,:) = {'lociXY', 'array', [0,0]};
defaults(end+1,:) = {'verbose','boolean',false};
defaults(end+1,:) = {'veryverbose','boolean',false};
defaults(end+1,:) = {'showplot','boolean',false};
defaults(end+1,:) = {'showPlots','boolean',false};
defaults(end+1,:) = {'showExtraPlots','boolean',false};
defaults(end+1,:) = {'overwrite','boolean',false};
defaults(end+1,:) = {'eTable','freeType',[]};

defaults(end+1,:) = {'saveData','boolean',false}; % changed default
defaults(end+1,:) = {'saveFolder','string',''};

% Parameters forwarded to ChrTracer3_FitSpots
% should use my parameter inheritances 
%   (I suppose that means these should really be classes that can actually
%   inherit from one another, rather than functions). 

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
            dataTable{s} = ChrTracer3_MultiFitSpots(fidSpts,datSpts,s,'parameters',pars);
        catch er
                warning(['encountered error on fov ',num2str(f), ' spot ',num2str(s)]);
                warning(er.getReport);
            if pars.stopOnError
                disp('set breakpoint here'); % note: matlab breakpoints don't work in parallel mode               
            end
        end
    end
else
    parfor s=selectSpots
        fidSpts = fidSpots{s};
        datSpts = datSpots{s};
        try
        dataTable{s} = ChrTracer3_MultiFitSpots(fidSpts,datSpts,s,'parameters',pars);
        catch 
           if pars.verbose
                warning(['encountered error on fov ',num2str(f), ' spot ',num2str(s)]);
                 % note: matlab breakpoints don't work in parallel mode               
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
