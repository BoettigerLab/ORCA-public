function datTable = ChrTracer3_FitSpots(fidSpts,datSpts,s,varargin)
% 
% fidSpts = cell array, Nhybs x 1, of 3D fiducial images
% datSpts = cell array, Nhybs x nDataChns, of 3D data images

% exTable =  E:\Alistair\2017-06-12_Emb10-12_BXC_cont\ExperimentLayout.xlsx
% saveFolder = E:\Alistair\2017-06-12_Emb10-12_BXC_cont_NewAnalysis\
% 
% 

defaults = cell(0,3);

% Fiducial alignment defaults
defaults(end+1,:) = {'refHybe','integer',1}; % hybe to use to start alignment
defaults(end+1,:) = {'upsample','positive',8};  % 8 for accuracy 2 for speed
defaults(end+1,:) = {'upsampleZ','positive',8};  % 8 for accuracy 2 for speed
defaults(end+1,:) = {'maxXYdrift','positive',4};
defaults(end+1,:) = {'maxZdrift','positive',6};
defaults(end+1,:) = {'fidRegTheta','nonnegative',.6};

% General FOV parameters
defaults(end+1,:) = {'fov','nonnegative',1};
defaults(end+1,:) = {'lociXY', 'array', [0,0]};
defaults(end+1,:) = {'saveFolder','string',''};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'veryverbose','boolean',false};
defaults(end+1,:) = {'showPlots','boolean',false};
defaults(end+1,:) = {'showExtraPlots','boolean',false};
defaults(end+1,:) = {'overwrite','boolean',false};
defaults(end+1,:) = {'eTable','freeType',[]};
defaults(end+1,:) = {'saveData','boolean',false}; % changed default

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
defaults(end+1,:) = {'datMinSep', 'nonnegative', 3};
defaults(end+1,:) = {'datMinHBratio','nonnegative',1.2}; % peak value over background value
defaults(end+1,:) = {'datMinAHratio','nonnegative',.25}; % fitted height over background vs peak value
defaults(end+1,:) = {'datMaxUncert','nonnegative',2}; % pixels
defaults(end+1,:) = {'keepBrightest','integer',1}; % will use variable idx to return the brightness rank of each spot  
defaults(end+1,:) = {'bkdFrac','fraction',0}; % 0 for off. dimmest bkdFrac of images will be used to compute local illumination background 


% parse defaults
pars = ParseVariableArguments(varargin,defaults,mfilename); 
%%
dispZpad = 2; % padding for fiducial z-display  (was 10)

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
    nChns = size(datSpts,2);
    [nRows,nCols,nStks] = size(datSpts{1}); % first hybe

    % Get some labels;
    if ~isempty(pars.eTable)
        datPropTable = DataChnsFromTable(pars.eTable);
    else % build it
        readout = (1:(numHybes*nChns))';
        chn = cellstr(repmat('data',numHybes*nChns,1));
        dataType = cellstr(repmat('H',numHybes*nChns,1));
        hybe = reshape(repmat(1:numHybes,nChns,1),numHybes*nChns,1);
        chnNum = reshape(repmat(1:nChns,1,numHybes),numHybes*nChns,1);
        datPropTable = table(readout,chn,dataType,hybe,chnNum); 
    end

    % Stack fiducials into a kernel. 
    % stackSpots = cat(4,fidSpts{:});
    % kernel = nanmedian(stackSpots,4);
    kernel = fidSpts{pars.refHybe}; % change behavior, use first spot. 
    % kernel = TranslateImage(kernel,1.3,1.3); % sanity check, manual offset

    % fit fiducial kernel to define restricted field of view
    if ~isinf(pars.fidMinPeakHeight) &&  pars.fidKeepBrightest > 0 
        fidTable = table();
        nTries = 5;
        n=0;
        while isempty(fidTable) && n<nTries   
            try
                fidTable = FindPeaks3D(kernel,...
                        'minPeakHeight',pars.fidMinPeakHeight,...
                        'maxFitWidth',pars.fidMaxFitWidth,...
                        'keepBrightest',1,... % pars.fidKeepBrightest
                        'minSep',pars.fidMinSep,...
                        'troubleshoot',pars.fidTroubleshoot); 
            catch
                fidTable = [];
                if pars.verbose
                    warning(['error processing spot ',num2str(s)]);
                    warning(['position ',num2str(lociXYx),' ',num2str(lociXYy)]);
                end
            end
            n=n+1;
            pars.fidMinPeakHeight = .6*pars.fidMinPeakHeight;
        end
        nFits = height(fidTable);
    
        % Define a highly restricted field of view to fit data in
        xs = cell(nFits,1);
        ys = cell(nFits,1);
        zs = cell(nFits,1);
        for i = 1:nFits % loop over all peaks in fiducial image
            % trim a box around the center of the fiducial image
            xi = max(round(fidTable.x(i)-pars.maxXYstep),1);
            xe = min(round(fidTable.x(i)+pars.maxXYstep),nCols);
            yi = max(round(fidTable.y(i)-pars.maxXYstep),1);
            ye = min(round(fidTable.y(i)+pars.maxXYstep),nRows);
            zi = max(round(fidTable.z(i)-pars.maxZstep), 1);
            ze = min(round(fidTable.z(i)+pars.maxZstep), nStks);
            xs{i} = xi:xe;
            ys{i} = yi:ye;
            zs{i} = zi:ze;
        end
        alignFid = true;
    else
        nFits = 1;
        xs{1} = 1:nCols;
        ys{1} = 1:nRows;
        zs{1} = 1:nStks;
        alignFid = false;
        x = 0; y=0; z=0; h=0; 
        fidTable = table(x,y,z,h);
    end
    %------------------ just plotting -----------------------------
    if pars.showPlots
        figure(2); clf;
        if alignFid
            stk = cellfun(@(x) max(x,[],3),fidSpts,'UniformOutput',false);
            stk0xy = cat(3,stk{:});
            stk = cellfun(@(x) max(permute(x,[3,2,1]),[],3),fidSpts,'UniformOutput',false);
            stk0xz = cat(3,stk{:});
            nZ = size(stk0xz,1); stk0xz = stk0xz(dispZpad:nZ-dispZpad,:,:);
            subplot(2,2,1);  Ncolor(2/numHybes*IncreaseContrast(stk0xy)); title('corr. fid x,y');
            hold on; plot(fidTable.x,fidTable.y,'r+');
            subplot(2,2,2);  Ncolor(2/numHybes*IncreaseContrast(stk0xz)); title('corr. fid x,z');
            hold on; plot(fidTable.x,fidTable.z-dispZpad,'r+');  
        end

        % Show absolute signal and background
        figure(10); clf;
        for d=1:nChns
            maxDat = cellfun(@(x) quantile(x(:),.999),datSpts(:,d));
            minDat = cellfun(@(x) quantile(x(:),.02),datSpts(:,d));
            subplot(nChns,1,d);
            bar(maxDat); hold on; bar(minDat);
            title(['channel ',num2str(d)]);
            legend('data spot','data background'); 
        end
        
        % BACKGROUND subtraction
        if pars.bkdFrac ~= 0
            for d=1:nChns % d = 2
                maxDat = cellfun(@(x) quantile(x(:),.999),datSpts(:,d));
                bkdHybIdx = (maxDat <= quantile(maxDat,pars.bkdFrac));
                bkdMap = nanmedian(cat(4,datSpts{bkdHybIdx,d}),4);
                figure(10+d); clf; ProjectIm3D(bkdMap);
                %datSptsOrig = datSpts;
                try  % handle errors on small datasets;
                    datSpts(:,d) = cellfun(@(x) x-bkdMap,datSpts(:,d),'UniformOutput',false);
                catch
                end        
            end
        end
    end
    %------------------------------------------------------------------------% 

    % Correct sub-pixel drift and fit data
    fidSptsAlign = cell(numHybes,1);
    fidSptsAlignTest = cell(numHybes,1);
    datSptsAlign = cell(numHybes,nChns);
    shifts = cell(numHybes,1);
    datTable = table();
    if ~isempty(fidTable.x)
        k=0;
        for h=1:numHybes % h=1
            % compute shifts (this should be accelerated)
            %   there is probably a more robust approach as well. 
            % fidSpts{1}
            if alignFid
                
                if pars.fidKeepBrightest == 1
                    ker = kernel;% ker = fidSpts{1}; % 
                    [fidSptsAlign{h},shifts{h}] = Register3D(ker,fidSpts{h},...
                        'center',[fidTable.x(1),fidTable.y(1),fidTable.z(1)],...
                        'threshold',pars.fidRegTheta,...
                        'upsample',pars.upsample,'verbose',pars.veryverbose,...
                        'maxShiftXY',pars.maxXYdrift,'maxShiftZ',pars.maxZdrift,...
                        'showplots',pars.showExtraPlots); % just for debugging
                elseif pars.fidKeepBrightest == 2
                     fidTable_i = FindPeaks3D(fidSpts{h},...
                        'minPeakHeight',pars.fidMinPeakHeight,...
                        'maxFitWidth',pars.fidMaxFitWidth,...
                        'keepBrightest',1,... % pars.fidKeepBrightest
                        'minSep',pars.fidMinSep,...
                        'troubleshoot',pars.fidTroubleshoot);
                    shifts{h}.xshift = -fidTable_i.x(1) + fidTable.x(1);
                    shifts{h}.yshift = -fidTable_i.y(1) + fidTable.y(1);
                    shifts{h}.zshift = -fidTable_i.z(1) + fidTable.z(1);
                    shifts{h}.score = 1;
                       if pars.veryverbose
                            disp(['computed ',...
                                ' xshift ', num2str(shifts{h}.xshift),...
                                ' yshift ',num2str(shifts{h}.yshift),...
                                ' zshift ',num2str(shifts{h}.zshift)]);
                       end
                elseif pars.fidKeepBrightest == 0
                    fidSptsAlign{h}=fidSpts{h};
                    shifts{h}.xshift = 0;
                    shifts{h}.yshift = 0;
                    shifts{h}.zshift = 0;
                    shifts{h}.score = 1;
                end
            
            % if pars.saveData || pars.showPlots
                % apply shifts to data channel
                for n=1:nChns
                    datSptsAlign{h,n} = TranslateImage( datSpts{h,n},shifts{h}.xshift,shifts{h}.yshift,...
                                    'zshift',shifts{h}.zshift,'upsample',pars.upsample,...
                                    'padValue',edgeValue(datSpts{h,n}));
                end
                fidSptsAlign{h}= TranslateImage( fidSpts{h},shifts{h}.xshift,shifts{h}.yshift,...
                                    'zshift',shifts{h}.zshift,'upsample',pars.upsample,...
                                    'padValue',edgeValue(fidSpts{h}));
            else
               datSptsAlign(h,:) = datSpts(h,:);
               shifts{h}.xshift = 0; 
               shifts{h}.yshift = 0;
               shifts{h}.zshift = 0;
               shifts{h}.score = 0;
            end
            % end
             % fit data channel in highly restricted field of view
            for n=1:nChns 
                k=k+1;
                for i=1:nFits 
                    try
                    % need to fit the aligned data if we are going to use
                    % the max distance from fiducial peak.
                    % could probably rewrite this to move the crop
                    % 
                    datSpot = datSptsAlign{h,n}(ys{i},xs{i},zs{i}); 
                    dTable = FitPsf3D(datSpot,...
                                'minPeakHeight',pars.datMinPeakHeight,...
                                'maxFitWidth',pars.datMaxFitWidth,...
                                'maxFitZdepth',pars.datMaxFitZdepth,...
                                'keepBrightest',pars.keepBrightest,...
                                'xyUnitConvert',pars.nmXYpix,...
                                'zUnitConvert',pars.nmZpix,...
                                'minHBratio',pars.datMinHBratio,...
                                'minAHratio',pars.datMinAHratio,...
                                'maxUncert',pars.datMaxUncert,...
                                'troubleshoot',pars.datTroubleshoot);   
                    % positions (add back the trim)
                    dTable.x = dTable.x + (xs{i}(1)-1)*pars.nmXYpix;  % in nm
                    dTable.y = dTable.y + (ys{i}(1)-1)*pars.nmXYpix;  % in nm
                    dTable.z = dTable.z + (zs{i}(1)-1)*pars.nmZpix;  % in nm
                    dTable.xL = dTable.xL + (xs{i}(1)-1)*pars.nmXYpix;  % in nm
                    dTable.yL = dTable.yL + (ys{i}(1)-1)*pars.nmXYpix;  % in nm
                    dTable.zL = dTable.zL + (zs{i}(1)-1)*pars.nmZpix;  % in nm
                    dTable.xU = dTable.xU + (xs{i}(1)-1)*pars.nmXYpix;  % in nm
                    dTable.yU = dTable.yU + (ys{i}(1)-1)*pars.nmXYpix;  % in nm
                    dTable.zU = dTable.zU + (zs{i}(1)-1)*pars.nmZpix;  % in nm
                    % from datPropTable
                    dTable.readout = datPropTable.readout(k)*ones(length(dTable.x),1);
                    dTable.dataType = repmat(string(datPropTable.dataType{k}),length(dTable.x),1); 
                    dTable.chn = repmat(string(datPropTable.chn{k}),length(dTable.x),1); 
                    % fov parameters
                    dTable.fov = pars.fov*ones(length(dTable.x),1);
                    dTable.locusX = lociXYx(s)*ones(length(dTable.x),1)*pars.nmXYpix;  % in nm
                    dTable.locusY = lociXYy(s)*ones(length(dTable.x),1)*pars.nmXYpix;  % in nm
                    % identifiers
                    dTable.tableOrder = h*ones(length(dTable.x),1);
                    dTable.hybe = datPropTable.hybe(k)*ones(length(dTable.x),1);
                    dTable.chnNum = n*ones(length(dTable.x),1);
                    dTable.idx = i*ones(length(dTable.x),1);
                    dTable.panel = (nChns*(h-1)+n)*ones(length(dTable.x),1);
                    dTable.s = s*ones(length(dTable.x),1);
                    dTable.fs = CantorPair(pars.fov,s)*ones(length(dTable.x),1);
                    dTable.fsr = CantorPair(dTable.fs,dTable.readout);
                    % fiducial corrections and fid positions
                    dTable.xshift = shifts{h}.xshift*ones(length(dTable.x),1)*pars.nmXYpix; % in nm
                    dTable.yshift = shifts{h}.yshift*ones(length(dTable.x),1)*pars.nmXYpix; % in nm
                    dTable.zshift = shifts{h}.zshift*ones(length(dTable.x),1)*pars.nmZpix; % in nm
                    dTable.shiftScore = shifts{h}.score*ones(length(dTable.x),1); % in nm
                    dTable.fid_x = fidTable.x(i)*ones(length(dTable.x),1)*pars.nmXYpix; % in nm
                    dTable.fid_y = fidTable.y(i)*ones(length(dTable.x),1)*pars.nmXYpix ;  % in nm 
                    dTable.fid_z = fidTable.z(i)*ones(length(dTable.x),1)*pars.nmZpix;  % in nm 
                    dTable.fid_h = fidTable.h(i)*ones(length(dTable.x),1);
%                     % corrected positions
%                     dTable.xc = dTable.x + dTable.xshift;
%                     dTable.yc = dTable.y + dTable.yshift;
%                     dTable.zc = dTable.z + dTable.zshift;
                    datTable = cat(1,datTable,dTable); 
                    catch er
                        if pars.verbose
                            warning(['encountered error on fov ',num2str(pars.fov), ' spot ',num2str(s), ' hyb ',num2str(h)]);
                            warning(er.getReport);
                        end
                    end
                end
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
        if pars.showPlots || pars.saveData
            % all this can be skipped for speed 
            
            tempDat = datSptsAlign';
            datMat = cat(4,tempDat{:});
            % datMat = cat(4,datSptsAlign{:});
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
                
                % get tile labels
                if ~isempty(pars.eTable)
                     [tileLabels_fid,tileLabels_dat] =TileLabelsFromEtable(pars.eTable);
                else
                    tileLabels_dat = cellstr(num2str( (1:numHybes*nChns)')) ;
                    tileLabels_fid = cellstr(num2str( (1:numHybes)')) ;
                end
                
                % 4D tiles
                if alignFid
                    figure(4); clf;
                    PlotProjection4D(fidMat,'fits',[],'projection','xy','tileLabels',tileLabels_fid); title('fid xy');
                    figure(5); clf;
                    PlotProjection4D(fidMat,'fits',[],'projection','xz','tileLabels',tileLabels_fid); title('fid xz');
                end
                figure(6); clf;
                PlotProjection4D(datMat,'fits',datTable,'projection','xy','tileLabels',tileLabels_dat,'nmXYpix',pars.nmXYpix,'nmZpix',pars.nmZpix); title('data xy');
                figure(7); clf;
                PlotProjection4D(datMat,'fits',datTable,'projection','xz','tileLabels',tileLabels_dat,'nmXYpix',pars.nmXYpix,'nmZpix',pars.nmZpix); title('data xz');

                % Overlay Projections
                if alignFid
                    figure(2);  % add corrected plots
                    stk = cellfun(@(x) max(x,[],3),fidSptsAlign,'UniformOutput',false);
                    stk1xy = cat(3,stk{:}); 
                    stk = cellfun(@(x) max(permute(x,[3,2,1]),[],3),fidSptsAlign,'UniformOutput',false);
                    stk1xz = cat(3,stk{:}); 
                    nZ = size(stk1xz,1); stk1xz = stk1xz(dispZpad:nZ-dispZpad,:,:);
                    subplot(2,2,3);  Ncolor(2/numHybes*IncreaseContrast(stk1xy)); title('corr. fid x,y aligned');
                    hold on; plot(fidTable.x,fidTable.y,'r+');
                    subplot(2,2,4); Ncolor(2/numHybes*IncreaseContrast(stk1xz)); title('corr. fid x,z  aligned');
                    hold on; plot(fidTable.x,fidTable.z-dispZpad,'r+');
                end

                % Data Projections
                isData = strcmp(datPropTable.dataType,'H');
                datSptsPlot = datSptsAlign(isData,:);
                figure(3);  clf;
                stk = cellfun(@(x) max(x,[],3),datSptsPlot,'UniformOutput',false);
                dat1xy = cat(3,stk{:});          
                stk = cellfun(@(x) max(permute(x,[3,2,1]),[],3),datSptsPlot,'UniformOutput',false);
                dat1xz = cat(3,stk{:});
                dat1xz(dat1xz==0) = mode(nonzeros(dat1xz(:)));
                subplot(1,2,1);  Ncolor(2/numHybes*IncreaseContrast(dat1xy)); title('corr. data x,y');
                subplot(1,2,2);  Ncolor(2/numHybes*IncreaseContrast(dat1xz)); title('corr. data x,z');
                % subplot(1,2,1); im = Ncolor(dat1xy); imagesc(IncreaseContrast(im)); title('corr. data x,y');
                % subplot(1,2,2); im = Ncolor(dat1xz); imagesc(IncreaseContrast(im)); title('corr. data x,z');
            end
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
    bx =  [i2(1:2,:,:), i2(end-1:end,:,:)];
    by =  [i2(:,1:2,:), i2(:,end-1:end,:)];
    bz =  [i2(:,:,1:2), i2(:,:,end-1:end)];
    edge = quantile( [bx(:) ; by(:); bz(:)],.8);