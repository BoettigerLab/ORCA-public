function spotTable = CT4_FitSpots(fidSptRef,fidSptCur,datSptCur,imProps,fidTable,varargin)
% fit spots

defaults = cell(0,3);
% Register3D: fine-scale drift correction with fiducial (sent to Register3D)
defaults(end+1,:) = {'threshold', 'fraction',0.6};
defaults(end+1,:) = {'upsample', 'integer',4};
defaults(end+1,:) = {'upsampleZ', 'integer',0};
defaults(end+1,:) = {'maxShiftXY', 'integer',8};
defaults(end+1,:) = {'maxShiftZ', 'integer',8};
defaults(end+1,:) = {'boxSize', 'integer',20};
defaults(end+1,:) = {'boxSizeZ', 'integer',20};
defaults(end+1,:) = {'showPlots', 'boolean',true};
defaults(end+1,:) = {'multiFit', 'boolean',true}; % switch to single fitter
defaults(end+1,:) = {'figShowAlign', 'freeType',12};
defaults(end+1,:) = {'verbose', 'boolean',true};

% FitOverlapSpots3D: data fitting
defaults(end+1,:) = {'nFits', 'integer',3};
defaults(end+1,:) = {'figShowFits', 'freeType',11};
% local 
% shared
defaults(end+1,:) = {'showExtraPlots', 'freeType',false};

pars = ParseVariableArguments(varargin,defaults,mfilename); 

if pars.upsampleZ == 0 % auto-assign to equivelant size
    xzRatio = imProps.xy2um/imProps.z2um;
    pars.upsampleZ = round(pars.upsample/xzRatio); 
end

% if pars.multiFit == false
%     pars.nFits = 1;
% end

currHyb = imProps.Hyb;
fovNum = imProps.fov;  % To check: this was failing to update properly...
nDataChns = length(datSptCur);
recordData = true; % useful for troubleshooting
    
    
    % --------- Initialize the spot table -------% 
    % this is the new spot table format. It is is a little bulky to write
    % it out like this, but it makes it obvious what all the variable names
    % are. Hard-coding it here is a conscious attempt to make the format
    % static and inflexible. 
    % ?? Maybe would have been a little more elegant to stack these in a
    %       structure,  st(N).spotID = [];  There are also no strings here,
    %       we just have ID numbers for barcodes that will need to be
    %       matched downstream. So matrices could be used, but they don't
    %       track the column heading as nicely.  
    N = pars.nFits*nDataChns;
    %=== Main table ================
    spotID = nan(N,1);
    traceID = nan(N,1);
    x_um = nan(N,1);  % spot position relative to fid (useful for distance filtering)  
    y_um = nan(N,1);  % spot position relative to fid
    z_um = nan(N,1);  % spot position relative to fid
    x_fid = nan(N,1); % fid position relative to FOV (useful for FOV specific errors like chromatic aberration and illumination)
    y_fid = nan(N,1); % fid position relative to FOV
    z_fid = nan(N,1); % fid position relative to FOV
    barcodeID = nan(N,1);
    cellID = nan(N,1);
    %= additional spot data==========
    fitQuality = nan(N,1); % false if fit outputs = inputs, denotes failed fit
    x_shift_um = nan(N,1);
    y_shift_um = nan(N,1);
    z_shift_um = nan(N,1);
    shift_score = nan(N,1);
    spt_brightness = nan(N,1);
    x_um_ci95 = nan(N,1);
    y_um_ci95 = nan(N,1);
    z_um_ci95 = nan(N,1);
    spt_bkd = nan(N,1);
    spt_bkd_xtilt = nan(N,1);
    spt_bkd_ytilt = nan(N,1);
    spt_sXY = nan(N,1);
    spt_sZ = nan(N,1);
    fid_brightness = nan(N,1);
    dat_chn = nan(N,1);
    % ---- data common to all spots
    stageX = nan(N,1); % fov position relative to stage/sample
    stageY = nan(N,1); % fov position relative to stage/sample
    fov  = nan(N,1);
    hyb = nan(N,1);
    % --------------------------------------------% 
    
     % Data associated with Ref Hyb does not need to be aligned
        % Z-crop
     if isempty(fidSptRef)
         fidSptRef = datSptCur{1};
     end
        zSteps = imProps.zSteps;
        [~,fidMaxIdx] = max(fidSptRef(:));
        [~,~,fidMax_z] = ind2sub(size(fidSptRef),fidMaxIdx);
        zi = max(1,fidMax_z-pars.boxSizeZ);
        ze = min(zSteps,fidMax_z+pars.boxSizeZ);
        % perform the registration
     if ~isempty(fidSptCur)
        [~,shifts] = Register3D( fidSptRef(:,:,zi:ze),...
                                 fidSptCur(:,:,zi:ze),...
                                'threshold',pars.fidRegTheta,...
                                'upsample',pars.upsample,...
                                'upsampleZ',pars.upsampleZ,...
                                'maxShiftXY',pars.maxXYdrift,...
                                'maxShiftZ',pars.maxZdrift,...
                                'showExtraPlots',false,...  %
                                'showplots',true,...  %
                                'figShowAlign',pars.figShowAlign,...
                                'verbose',pars.verbose); % just for debugging
    %    %  just for troubleshooting
    %     test = imtranslate(fidSptCur(:,:,zi:ze),[shifts.xshift,shifts.yshift,shifts.zshift]);
    %     figure(10); clf;
    %     [test_xy,test_xz] = ProjectIm3D(test,'showPlots',true);
    %     [orig_xy,orig_xz] = ProjectIm3D(fidSptRef(:,:,zi:ze),'showPlots',true);
    %     figure(10); clf; subplot(1,2,1); Ncolor(IncreaseContrast(cat(3,orig_xy,test_xy)));
    %      subplot(1,2,2); Ncolor(IncreaseContrast(cat(3,orig_xz,test_xz)));
    else
       shifts.xshift=0; 
       shifts.yshift=0;
       shifts.zshift=0;
       shifts.score =1;
    end

    %% Fitting
    %   For each trace, s,
    %   first take the n-brightest pixels in each barcode 
    %   then go back and fit 3D PSFs centered on each of these
    %   record all n for each hyb, even if we have to populate with nans, this
    %        will make the indexing easier
    %   then take the brightest spot in each barcode to define a 'cloud' and a
    %       centroid for the trace. 
    %   for each of the 5 spots, we have now 3 quality values: brightness,
    %       distance-from-centroid, quality-of-PSF-fit, per trace 
    %   these quality values can be compared across all traces for the same
    %       hyb, to select a maximum likelihood candidate among them. One could
    %       use 'brightest' fraction to create an estimated distribution over
    %       distances, or a                  
    bx = floor(pars.maxXYdrift/2);
    bz = floor(pars.maxZdrift/2);
    
    for d=1:nDataChns
        xErrWidth = nan; % place holder
        currSpotIm = datSptCur{d}(:,:,zi:ze);
        bkd = median(currSpotIm(:));
        currSpotIm = imtranslate(currSpotIm,[shifts.xshift,shifts.yshift,shifts.zshift],'FillValues',bkd);
        data_3d = currSpotIm(bx:end-bx,bx:end-bx,bz:end-bz);

        % seed with brightest non-joined pixels
        bw = imregionalmax(data_3d); % 26-connectivity is default for 3D
        dataLocalMax = data_3d;
        dataLocalMax(~bw) = 0;
        nFits = pars.nFits;
        [~,id] = sort(dataLocalMax(:),'descend') ; % data_3d = im3D;
        x = zeros(nFits,1);
        y = zeros(nFits,1);
        z = zeros(nFits,1);
        h = zeros(nFits,1);
        for n=1:nFits
            [y(n),x(n),z(n)] = ind2sub(size(data_3d),id(n));
            h(n) = data_3d(y(n),x(n),z(n));
        end   
        if pars.showExtraPlots
            figure(10); clf; 
            ProjectIm3D(data_3d); colormap(gray);
            subplot(1,3,1); hold on; plot(x,y,'r.');
            subplot(1,3,2); hold on; plot(y,z,'r.');
            subplot(1,3,3); hold on; plot(x,z,'r.');
            pause(.1);
        end
        %%
        % do we fit these as N overlapping PSFs? or % 
        %   3 is much faster than 5. 
        %   if we don't fit as overlapping spots the single bright image is going
        %   to drag the others in with it.  
        %   (Also a question to figure out with Bogdan's version).
        % adjust expected PSF sizes based on voxel size
        minSigma = 0.9*.1/imProps.xy2um;
        maxSigma = 2.4*.1/imProps.xy2um;
        initSigma= 1.5*.1/imProps.xy2um;
        minSigmaZ = 0.9*.1/imProps.z2um*3;
        maxSigmaZ = 2.4*.1/imProps.z2um*3;
        initSigmaZ= 1.5*.1/imProps.z2um*3;
        if pars.multiFit  % also just pass it nFits =1
        [fitTable,fitError] = FitOverlapSpots3D(data_3d,nFits,'initXYZ',[x,y,z],...
            'showplot',pars.showPlots,'figShowFits',pars.figShowFits,...
            'brightVar',3,'maxTilt',100,...
            'minSigma',minSigma,'maxSigma',maxSigma,'initSigma',initSigma,...
            'minSigmaZ',minSigmaZ,'maxSigmaZ',maxSigmaZ,'initSigmaZ',initSigmaZ);
        else
          fitError = [];
          fitTable = FitPsf3D(data_3d,...
                                'seedPoint',[x,y,z,h],...
                                'minPeakHeight',200,...
                                'maxFitWidth',8,...
                                'maxFitZdepth',8,...
                                'keepBrightest',nFits,...
                                'xyUnitConvert',1,...  108  % converted below
                                'zUnitConvert',1,... 250
                                'minHBratio',.6,...
                                'minAHratio',.2,...
                                'maxUncert',2,...
                                'maxSigma',1.5,...
                                'initSigmaXY',.7,...
                                'troubleshoot',pars.showPlots,...
                                'filterSpot',false);    % the whole point is to filter later
                            % convert some names
                            fitTable.sXY = fitTable.wx;
                            fitTable.sZ = fitTable.wz;
                            fitTable.px = zeros(size(fitTable.x));
                            fitTable.py = zeros(size(fitTable.x));
                            xErrWidth = fitTable.xU - fitTable.xL;
                            yErrWidth = fitTable.yU - fitTable.yL;
                            zErrWidth = fitTable.zU - fitTable.zL;
                            % a and b are defined the same
            
        end
        %% ============= organize into table ====================
        if ~isempty(fitError) && isnan(xErrWidth)
            xErrWidth = fitError.ci95(6:4:end,2)-fitError.ci95(6:4:end,1);
            yErrWidth = fitError.ci95(7:4:end,2)-fitError.ci95(7:4:end,1);
            zErrWidth = fitError.ci95(8:4:end,2)-fitError.ci95(8:4:end,1);
        end 
        xy2um = imProps.xy2um;
        z2um = imProps.z2um;
        
        if recordData
            % spotID traceID x_um y_um z_um barcodeID cellID (barcodeID will become chr,start,end)              
            for nn=1:nFits
                n = (d-1)*nFits+nn;
                % this is the main data table
                spotID(n) = nan; % will be populated in master table;
                traceID(n) = fidTable.traceID; %  
                x_um(n) = (-pars.boxSize+bx + fitTable.x(nn))*xy2um;
                y_um(n) = (-pars.boxSize+bx + fitTable.y(nn))*xy2um;
                z_um(n) = (fitTable.z(nn)+bz)*z2um;
                barcodeID(n) = nDataChns*(currHyb-1)+d;  % load from experiment table
                cellID(n) = fidTable.cellID;
                % this is the additional data table  
                fitQuality(n) = fitTable.fitQuality(nn);
                x_fid(n) = fidTable.x_pix*xy2um;
                y_fid(n) = fidTable.y_pix*xy2um;
                z_fid(n) = (zi+pars.boxSizeZ)*z2um;
                x_shift_um(n) = shifts.xshift*xy2um;
                y_shift_um(n) = shifts.yshift*xy2um;
                z_shift_um(n) = shifts.zshift*z2um;
                shift_score(n) = shifts.score;
                spt_brightness(n) = fitTable.a(nn);
                x_um_ci95(n) = xErrWidth(nn)*xy2um;
                y_um_ci95(n) = yErrWidth(nn)*xy2um;
                z_um_ci95(n) = zErrWidth(nn)*z2um;
                spt_bkd(n) = fitTable.b(nn);
                spt_bkd_xtilt(n) = fitTable.px(nn);
                spt_bkd_ytilt(n) = fitTable.py(nn);
                spt_sXY(n) = fitTable.sXY(nn);
                spt_sZ(n) = fitTable.sZ(nn);
                fid_brightness(n) = fidTable.brightness;
                dat_chn(n) = d;
                % constant data
                % these are actually constant for the same raw image-stack (same fov + hyb)  
                % it is convienent to have these in the same table though. 
                stageX(n) = imProps.stageXY(1); % different stages may report different units (mm, um)? 
                stageY(n) = imProps.stageXY(2);
                fov(n) = fidTable.fov;%  f; 
                hyb(n) = currHyb;
            end
        end
    end % end loop over data channels
spotTable = table(spotID,traceID,x_um,y_um,z_um,barcodeID,cellID,...
        fitQuality,x_fid,y_fid,z_fid,x_shift_um,y_shift_um,z_shift_um,shift_score,spt_brightness,x_um_ci95,y_um_ci95,z_um_ci95,...
        spt_bkd,spt_bkd_xtilt,spt_bkd_ytilt,spt_sXY,spt_sZ,fid_brightness,dat_chn,...
        stageX,stageY,fov,hyb);
    % end loop over spots in FOV  
    % end loop over FOVs 

    