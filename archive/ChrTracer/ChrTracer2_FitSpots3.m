function datTable = ChrTracer2_FitSpots3(fidSpts,datSpts,s,varargin)
% Align using max-projection of xy and xz fiducial images rather than the
% 3D correlation function used in ChrTracer2_FitSpots 
%
% exTable =  E:\Alistair\2017-06-12_Emb10-12_BXC_cont\ExperimentLayout.xlsx
% saveFolder = E:\Alistair\2017-06-12_Emb10-12_BXC_cont_NewAnalysis\
% ChrTracer2

defaults = cell(0,3);

% Fiducial alignment defaults
defaults(end+1,:) = {'upsample','positive',4};  % 8 for accuracy 2 for speed
defaults(end+1,:) = {'maxShiftXY','positive',4};
defaults(end+1,:) = {'maxShiftZ','positive',6};
defaults(end+1,:) = {'fidRegTheta','nonnegative',.6};

% General FOV parameters
defaults(end+1,:) = {'nmXYpix','positive',154};
defaults(end+1,:) = {'nmZpix','positive',100};
defaults(end+1,:) = {'fov','nonnegative',1};
defaults(end+1,:) = {'lociXY', 'array', [0,0]};
defaults(end+1,:) = {'saveFolder','string',''};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'veryverbose','boolean',false};
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
% pars = ParseVariableArguments([],defaults,[]);
pars = ParseVariableArguments(varargin,defaults,mfilename); 
%%
% global CT
% loading test data straight from CT
% fidSpots = CT{1}.fidSpots;
% datSpots = CT{1}.dataSpots;
% pars.lociXY = CT{1}.parsFOV.lociXY;
% pars.fov = CT{1}.parsFOV.fov;
% pars.dataFolder = CT{1}.parsFOV.saveFolder;


%% Setup file save names
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

if ~(exist(tableSaveName,'file')==2 && exist(imageSaveName,'file')==2) || pars.overwrite

    numHybes = length(fidSpts);
    [nRows,nCols,nStks] = size(datSpts{1}); % first hybe

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

%     % just for testing
%     fidSpts = {};
%     datSpts = {};
%     for i=1:size(fidMat,4)
%         fidSpts{i} = fidMat(:,:,:,i);
%         datSpts{i} = datMat(:,:,:,i);
%     end
%     
    % Stack fiducials into a kernel. 
    stackSpots = cat(4,fidSpts{:});
    kernel = nanmedian(stackSpots,4);
    % kernel = TranslateImage(kernel,1.3,1.3); % sanity check, manual offset

    % fit fiducial kernel to define restricted field of view
    fidTable = table();
    nTries = 5;
    n=0;
    while isempty(fidTable) && n<nTries      
        fidTable = FindPeaks3D(kernel,...
                'minPeakHeight',pars.fidMinPeakHeight,...
                'maxFitWidth',pars.fidMaxFitWidth,...
                'keepBrightest',pars.fidKeepBrightest,...
                'minSep',pars.fidMinSep,...
                'troubleshoot',pars.fidTroubleshoot); 
        n=n+1;
        pars.fidMinPeakHeight = .6*pars.fidMinPeakHeight;
    end
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

    %------------------ just plotting -----------------------------
    if pars.showPlots
        figure(2); clf;
        stk = cellfun(@(x) max(x,[],3),fidSpts,'UniformOutput',false);
        stk0xy = cat(3,stk{:});
        stk = cellfun(@(x) max(permute(x,[3,2,1]),[],3),fidSpts,'UniformOutput',false);
        stk0xz = cat(3,stk{:});
        subplot(2,2,1); im = Ncolor(2/numHybes*stk0xy); imagesc(IncreaseContrast(im)); title('corr. fid x,y');
        hold on; plot(fidTable.x,fidTable.y,'r+');
        subplot(2,2,2); im = Ncolor(2/numHybes*stk0xz); imagesc(IncreaseContrast(im)); title('corr. fid x,z');
        hold on; plot(fidTable.x,fidTable.z,'r+');
        
    end
    %------------------------------------------------------------------------% 

    % Correct sub-pixel drift and fit data
    fidSptsAlign = cell(numHybes,1);
    fidSptsAlignTest = cell(numHybes,1);
    datSptsAlign = cell(numHybes,1);
    shifts = cell(numHybes,1);
    datTable = table();
    if ~isempty(fidTable.x)
        for h=1:numHybes % h=1
            
            refXY = imresize( max(fidSpts{pars.refHyb},[],3), pars.upsample);  % could have a flag to use average spot as ref
            inpXY = imresize( max(fidSpts{h},[],3), pars.upsample);
            [xshift,yshift,angle] = CorrAlignRotate(refXY,inpXY,'angles',0)
            fidSptsAlign{h} = TranslateImage(inpXY,xshift,yshift);
        end
        
        stk = cellfun(@(x) max(x,[],3),fidSptsAlign,'UniformOutput',false);
        stk0xy = cat(3,stk{:});
        figure(3); clf; subplot(2,2,1); im = Ncolor(2/numHybes*stk0xy); imagesc(IncreaseContrast(im)); title('corr. fid x,y');
        
        
        for h=1:numHybes % h=1
            
            refXY = imresize( max(fidSpts{pars.refHyb},[],3), pars.upsample);  % could have a flag to use average spot as ref
            inpXY = imresize( max(fidSpts{h},[],3), pars.upsample);
            [xshift,yshift,angle] = CorrAlignRotate(refXY,inpXY,'angles',0);
            fidSptsAlign{h} = TranslateImage(inpXY,xshift,yshift);
        end
            

            % apply shifts to data channel
            datSptsAlign{h} = TranslateImage(datSpts{h},shifts{h}.xshift,shifts{h}.yshift,...
                                'zshift',shifts{h}.zshift,'upsample',pars.upsample,...
                                'padValue',edgeValue(datSpts{h}));
            fidSptsAlignTest{h}= TranslateImage(fidSpts{h},shifts{h}.xshift,shifts{h}.yshift,...
                                'zshift',shifts{h}.zshift,'upsample',pars.upsample,...
                                'padValue',edgeValue(fidSpts{h}));
             % fit data channel in highly restricted field of view
            for i=1:nFits
                datSpot = datSptsAlign{h}(ys{i},xs{i},zs{i}); 
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
                dTable.x = dTable.x + (xs{i}(1)-1)*pars.nmXYpix;  % in nm
                dTable.y = dTable.y + (ys{i}(1)-1)*pars.nmXYpix;  % in nm
                dTable.z = dTable.z + (zs{i}(1)-1)*pars.nmZpix;  % in nm
                dTable.read = bitNum(h)*ones(length(dTable.x),1);
                dTable.dataType = repmat(string(dataType{h}),length(dTable.x),1); 
                dTable.hybe = h*ones(length(dTable.x),1);
                dTable.idx = i*ones(length(dTable.x),1);
                dTable.fov = pars.fov*ones(length(dTable.x),1);
                dTable.locusX = lociXYx(s)*ones(length(dTable.x),1);
                dTable.locusY = lociXYy(s)*ones(length(dTable.x),1);
                dTable.s = s*ones(length(dTable.x),1);
                dTable.xshift = shifts{h}.xshift*ones(length(dTable.x),1); % in pix
                dTable.yshift = shifts{h}.yshift*ones(length(dTable.x),1); % in pix
                dTable.zshift = shifts{h}.zshift*ones(length(dTable.x),1); % in pix
                dTable.fid_x = fidTable.x(i)*ones(length(dTable.x),1)*pars.nmXYpix; % in nm
                dTable.fid_y = fidTable.y(i)*ones(length(dTable.x),1)*pars.nmXYpix ;  % in nm 
                dTable.fid_z = fidTable.z(i)*ones(length(dTable.x),1)*pars.nmZpix;  % in nm 
                dTable.fid_h = fidTable.h(i)*ones(length(dTable.x),1);
                datTable = cat(1,datTable,dTable); 
            end
        end
        %% save data table to disk

        if pars.saveData
            writetable(datTable,tableSaveName); 
            if pars.verbose
                disp(['wrote ',tableSaveName]);
            end
        end

        %% save data in binary format

        datMat = cat(4,datSptsAlign{:});
        description = ['channel=data',...
            '; fov=',num2str(pars.fov),'; spot=',num2str(s),...
            '; locus=(',num2str(lociXYx(s)),',',num2str(lociXYy(s)),')',...
            '; dataFolder=',pars.saveFolder];
        if pars.saveData
            SaveImage4D(datMat,imageSaveName,'description',description);
        end

        fidMat = cat(4,fidSptsAlign{:});
        description = ['channel=fiducial',...
            '; fov=',num2str(pars.fov),'; spot=',num2str(s),...
            '; locus=(',num2str(lociXYx(s)),',',num2str(lociXYy(s)),')',...
            '; dataFolder=',pars.saveFolder];
        if pars.saveData
            SaveImage4D(fidMat,fidSaveName,'description',description);
        end


        %% for troubleshooting
        if pars.showPlots 
            figure(4); clf;
            PlotProjection4D(fidMat,'fits',[],'projection','xy');
            figure(5); clf;
            PlotProjection4D(fidMat,'fits',[],'projection','xz');

            figure(6); clf;
            PlotProjection4D(datMat,'fits',datTable,'projection','xy');
            figure(7); clf;
            PlotProjection4D(datMat,'fits',datTable,'projection','xz');

            % Projections
            figure(2);  % add corrected plots
            stk = cellfun(@(x) max(x,[],3),fidSptsAlign,'UniformOutput',false);
            stk1xy = cat(3,stk{:});
            stk = cellfun(@(x) max(permute(x,[3,2,1]),[],3),fidSptsAlign,'UniformOutput',false);
            stk1xz = cat(3,stk{:});
            subplot(2,2,3); im = Ncolor(stk1xy); imagesc(IncreaseContrast(im)); title('corr. fid x,y');
            hold on; plot(fidTable.x,fidTable.y,'r+');
            subplot(2,2,4); im = Ncolor(stk1xz); imagesc(IncreaseContrast(im)); title('corr. fid x,z');
            hold on; plot(fidTable.x,fidTable.z,'r+');

            % Data Projections
            figure(3);  clf;
            stk = cellfun(@(x) max(x,[],3),datSptsAlign,'UniformOutput',false);
            dat1xy = cat(3,stk{:});
            stk = cellfun(@(x) max(permute(x,[3,2,1]),[],3),datSptsAlign,'UniformOutput',false);
            dat1xz = cat(3,stk{:});
            subplot(1,2,1); im = Ncolor(dat1xy); imagesc(IncreaseContrast(im)); title('corr. data x,y');
            subplot(1,2,2); im = Ncolor(dat1xz); imagesc(IncreaseContrast(im)); title('corr. data x,z');
        end

    else
        if pars.verbose
           disp('fiducial not above min intensity'); 
        end
    end
else
    if pars.verbose
       disp(['found existing data for spot ',num2str(s), '. Skipping to next spot.']);
    end
    datTable = readtable(tableSaveName);
end


function edge = edgeValue(i2)
    bx =  [i2(1,:,:), i2(end,:,:)];
    by =  [i2(:,1,:), i2(:,end,:)];
    bz =  [i2(:,:,1), i2(:,:,end)];
    edge = quantile( [bx(:) ; by(:); bz(:)],.8);