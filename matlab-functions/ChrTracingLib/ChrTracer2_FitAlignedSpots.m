function datTable = ChrTracer2_FitAlignedSpots(fidMat,datMat,s,varargin)

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
defaults(end+1,:) = {'saveData','boolean',true};

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
defaults(end+1,:) = {'tag','string','NewFits'}; % pixels
defaults(end+1,:) = {'shifts','freeType',[]}; % pixels
defaults(end+1,:) = {'fidGain','positive',.75};
defaults(end+1,:) = {'datGain','positive',.75};
defaults(end+1,:) = {'saveTable','boolean',true};

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


if sum(pars.lociXY) == 0 
    lociXYx = nan;
    lociXYy = nan;
else   
    lociXYx = pars.lociXY(s,1);
    lociXYy = pars.lociXY(s,2);
end
tableSaveName = [pars.saveFolder, 'fov',num2str(pars.fov,'%03d'),...
    '_spot',num2str(s,'%04d'),...
    '_locus(',num2str(lociXYx),',',num2str(lociXYy),')','_',pars.tag,'.csv'];

if ~(exist(tableSaveName,'file')==2)  || pars.overwrite

    
    [nRows,nCols,nStks,numHybes] = size(fidMat); % first hybe
    % Get some labels;
    if ~isempty(pars.eTable)
        dataType = pars.eTable.DataType;
        if ~iscell(pars.eTable.Bit(1))
            bitNum = pars.eTable.Bit;
        else
            bitNum = cellfun(@(x) str2double(strsplit(x,',')),pars.eTable.Bit,'UniformOutput',false); % added st2double
            bitNum = cat(1,bitNum{:});
        end
    else
        bitNum = 1:numHybes;
    end
    if isempty(pars.shifts)
       shifts = pars.shifts;
    else
        shifts = [];
    end

    % Stack fiducials into a kernel. 
    kernel = nanmedian(fidMat,4);
    % kernel = TranslateImage(kernel,1.3,1.3); % sanity check, manual offset

    % fit fiducial kernel to define restricted field of view
    fidTable = FindPeaks3D(kernel,...
                'minPeakHeight',pars.fidMinPeakHeight,...
                'maxFitWidth',pars.fidMaxFitWidth,...
                'keepBrightest',pars.fidKeepBrightest,...
                'troubleshoot',pars.fidTroubleshoot);  
    nFits = height(fidTable);
    % Define a highly restricted field of view to fit data in
    xs = cell(nFits,1);
    ys = cell(nFits,1);
    zs = cell(nFits,1);
    for i = 1:nFits % loop over all peaks in fiducial image
        xi = max(round(fidTable.x(i)-pars.datBoxXY),1);
        xe = min(round(fidTable.x(i)+pars.datBoxXY),nCols);
        yi = max(round(fidTable.y(i)-pars.datBoxXY),1);
        ye = min(round(fidTable.y(i)+pars.datBoxXY),nRows);
        zi = max(round(fidTable.z(i)-pars.datBoxZ), 1);
        ze = min(round(fidTable.z(i)+pars.datBoxZ), nStks);
        xs{i} = xi:xe;
        ys{i} = yi:ye;
        zs{i} = zi:ze;
    end

    % Correct sub-pixel drift and fit data
    datTable = table();
    for h=1:numHybes % h=1
         % fit data channel in highly restricted field of view
        for i=1:nFits
            datSpot = datMat(ys{i},xs{i},zs{i},h); 
            dTable = FitPsf3D(datSpot,...
                        'minPeakHeight',pars.datMinPeakHeight,...
                        'maxFitWidth',pars.datMaxFitWidth,...
                        'keepBrightest',1,...
                        'xyUnitConvert',pars.nmXYpix,...
                        'zUnitConvert',pars.nmZpix,...
                        'minHBratio',pars.datMinHBratio,...
                        'minAHratio',pars.datMinAHratio,...
                        'maxUncert',pars.datMaxUncert,...
                        'troubleshoot',pars.datTroubleshoot);         
            dTable.x = dTable.x + (xs{i}(1)-1)*pars.nmXYpix;
            dTable.y = dTable.y + (ys{i}(1)-1)*pars.nmXYpix;
            dTable.z = dTable.z + (zs{i}(1)-1)*pars.nmZpix;
            dTable.read = bitNum(h)*ones(length(dTable.x),1);
            dTable.dataType = repmat(string(dataType{h}),length(dTable.x),1); 
            dTable.hybe = h*ones(length(dTable.x),1);
            dTable.idx = i*ones(length(dTable.x),1);
            dTable.fov = pars.fov*ones(length(dTable.x),1);
            dTable.locusX = lociXYx*ones(length(dTable.x),1);
            dTable.locusY = lociXYy*ones(length(dTable.x),1);
            dTable.s = s*ones(length(dTable.x),1);
            dTable.fid_x = fidTable.x(i)*ones(length(dTable.x),1);
            dTable.fid_y = fidTable.y(i)*ones(length(dTable.x),1);
            dTable.fid_z = fidTable.z(i)*ones(length(dTable.x),1);
            dTable.fid_h = fidTable.h(i)*ones(length(dTable.x),1);
            if ~isempty(shifts)
                dTable.xshift = shifts{h}.xshift*ones(length(dTable.x),1); 
                dTable.yshift = shifts{h}.yshift*ones(length(dTable.x),1);
                dTable.zshift = shifts{h}.zshift*ones(length(dTable.x),1);
            end
            datTable = cat(1,datTable,dTable); 
        end
    end
    %% save data table to disk

    if pars.saveData && pars.saveTable
        writetable(datTable,tableSaveName); 
        if pars.verbose
            disp(['wrote ',tableSaveName]);
        end
    end

    %% for troubleshooting
    if isempty(datTable)
        if pars.verbose
            disp(['no fits recovered for spot: ',num2str(s)]);
        end
    elseif pars.showPlots  
        figure(1); clf;
        is1 = datTable.idx == 1;
        is2 = datTable.idx == 2;
        fidXY = squeeze(max(fidMat,[],3));
        fidXZ = squeeze(max(permute(fidMat,[3,2,1,4]),[],3));
        datXY = squeeze(max(datMat,[],3));
        datXZ = squeeze(max(permute(datMat,[3,2,1,4]),[],3));   
        subplot(2,2,1); Ncolor(pars.fidGain*fidXY);
            hold on; plot(fidTable.x, fidTable.y,'yo');
        subplot(2,2,2); Ncolor(pars.datGain*datXY);
            hold on; scatter(datTable.x(is1)/pars.nmXYpix,datTable.y(is1)/pars.nmXYpix,12,hsv(length(datTable.x(is1))),'o'  );
            hold on; scatter(datTable.x(is2)/pars.nmXYpix,datTable.y(is2)/pars.nmXYpix,20,hsv(length(datTable.x(is2))) ,'s' );
        subplot(2,2,3); Ncolor(pars.fidGain*fidXZ);
            hold on; plot(fidTable.x, fidTable.z,'yo');
        subplot(2,2,4); Ncolor(pars.datGain*datXZ);
            hold on; scatter(datTable.x(is1)/pars.nmXYpix,datTable.z(is1)/pars.nmZpix,12,hsv(length(datTable.x(is1))),'o'  );
            hold on; scatter(datTable.x(is2)/pars.nmXYpix,datTable.z(is2)/pars.nmZpix,20,hsv(length(datTable.x(is2))) ,'s' );
        
        figure(2); clf;
        PlotProjection4D(fidMat,'fits',datTable,'projection','xy');
        figure(3); clf;
        PlotProjection4D(fidMat,'fits',datTable,'projection','xz');
        
        figure(4); clf;
        PlotProjection4D(datMat,'fits',datTable,'projection','xy');
        figure(5); clf;
        PlotProjection4D(datMat,'fits',datTable,'projection','xz');  
    end
else
    if pars.verbose
       disp(['found existing data for spot ',num2str(s), '. Skipping to next spot.']);
    end
    datTable = readtable(tableSaveName);
end