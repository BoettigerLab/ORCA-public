function datTable = ChrTracer2_FitSpots(fidSpts,datSpts,s,varargin)

% exTable =  E:\Alistair\2017-06-12_Emb10-12_BXC_cont\ExperimentLayout.xlsx
% saveFolder = E:\Alistair\2017-06-12_Emb10-12_BXC_cont_NewAnalysis\
% ChrTracer2

defaults = cell(0,3);

% Fiducial alignment defaults
defaults(end+1,:) = {'upsample','positive',4};  % 8 for accuracy 2 for speed
defaults(end+1,:) = {'maxShiftXY','positive',4};
defaults(end+1,:) = {'maxShiftZ','positive',6};

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

%-----------(These just get forwarded to the fitting program)------------%
% Fiducial fitting defaults
defaults(end+1,:) = {'fidMinPeakHeight', 'positive', 200};
defaults(end+1,:) = {'fidCameraBackground', 'nonnegative', 0};
defaults(end+1,:) = {'fidPeakBlur', 'nonnegative', .5};
defaults(end+1,:) = {'fidTroubleshoot', 'boolean', false};
defaults(end+1,:) = {'fidMaxFitWidth', 'positive', 8};
defaults(end+1,:) = {'fidMinSep', 'nonnegative', 5};  % Min separation between peaks in pixels.  Closer than this will be averaged
defaults(end+1,:) = {'fidKeepBrightest','integer',2};  % Max number of peaks to allow
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
% Fit Aligned Spots specific
defaults(end+1,:) = {'tag','string','fits'}; % pixels

% parse defaults
pars = ParseVariableArguments(varargin,defaults,mfilename); 
%%
% global CT
% loading test data straight from CT
% fidSpots = CT{1}.fidSpots;
% datSpots = CT{1}.dataSpots;
% pars.lociXY = CT{1}.parsFOV.lociXY;
% pars.fov = CT{1}.parsFOV.fov;
% pars.dataFolder = CT{1}.parsFOV.saveFolder;

%  exist(tableSaveName,'file')==2

lociXYx = pars.lociXY(:,1);
lociXYy = pars.lociXY(:,2);
tableSaveName = [pars.saveFolder, 'fov',num2str(pars.fov,'%03d'),...
    '_spot',num2str(s,'%04d'),...
    '_locus(',num2str(lociXYx(s)),',',num2str(lociXYy(s)),')','_fits.csv'];
imageSaveName = [pars.saveFolder, 'fov',num2str(pars.fov,'%03d'),...
    '_spot',num2str(s,'%04d'),...
    '_locus(',num2str(lociXYx(s)),',',num2str(lociXYy(s)),')','_AlignedData.i4d'];
fidSaveName = [pars.saveFolder, 'fov',num2str(pars.fov,'%03d'),...
    '_spot',num2str(s,'%04d'),...
    '_locus(',num2str(lociXYx(s)),',',num2str(lociXYy(s)),')','_AlignedFid.i4d'];

if ~( exist(fidSaveName,'file')==2   && exist(imageSaveName,'file')==2   ) || pars.overwrite

    numHybes = length(fidSpts);   
    
    % Stack fiducials into a kernel. 
    stackSpots = cat(4,fidSpts{:});
    kernel = nanmedian(stackSpots,4);
    % kernel = TranslateImage(kernel,1.3,1.3); % sanity check, manual offset

    % Correct sub-pixel drift and fit data
    fidSptsAlign = cell(numHybes,1);
    fidSptsAlignTest = cell(numHybes,1);
    datSptsAlign = cell(numHybes,1);
    shifts = cell(numHybes,1);
    datTable = table();
    for h=1:numHybes % h=1
        % compute shifts (this should be accelerated)
        [fidSptsAlign{h},shifts{h}] = Register3D(kernel,fidSpts{h},...
            'upsample',pars.upsample,'verbose',false,...
            'maxShiftXY',pars.maxShiftXY,'maxShiftZ',pars.maxShiftZ);
        % apply shifts to data channel
        datSptsAlign{h} = TranslateImage(datSpts{h},shifts{h}.xshift,shifts{h}.yshift,...
                            'zshift',shifts{h}.zshift,'upsample',pars.upsample);
        fidSptsAlignTest{h}= TranslateImage(fidSpts{h},shifts{h}.xshift,shifts{h}.yshift,...
                            'zshift',shifts{h}.zshift,'upsample',pars.upsample);
    end

    %% save data in binary format
    datMat = cat(4,datSptsAlign{:});
    description = ['channel=data',...
        '; fov=',num2str(pars.fov),'; spot=',num2str(s),...
        '; locus=(',num2str(lociXYx(s)),',',num2str(lociXYy(s)),')',...
        '; dataFolder=',pars.saveFolder];
    SaveImage4D(datMat,imageSaveName,'description',description);

    fidMat = cat(4,fidSptsAlign{:});
    description = ['channel=fiducial',...
        '; fov=',num2str(pars.fov),'; spot=',num2str(s),...
        '; locus=(',num2str(lociXYx(s)),',',num2str(lociXYy(s)),')',...
        '; dataFolder=',pars.saveFolder];
    SaveImage4D(fidMat,fidSaveName,'description',description);

    
    %% Now we fit spots
    % this is a separate functions so we can call it later.
    datTable = ChrTracer2_FitAlignedSpots(fidMat,datMat,s,'parameters',pars,'shifts',shifts);
    
else
    if pars.verbose
       disp(['found existing data for spot ',num2str(s), '. Skipping to next spot.']);
    end
    datTable = readtable(tableSaveName);
end