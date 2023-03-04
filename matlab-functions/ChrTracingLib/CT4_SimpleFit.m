function spotTable = CT4_SimpleFit(cropIm,imProps,sptPixTable,d,nDataChns,currHyb,varargin)
% a simple wrapper of FitPsf3D, which arranges the data into the new table
% format.  

defaults = cell(0,3);
defaults(end+1,:) = {'boxSize', 'integer',4};
defaults(end+1,:) = {'boxSizeZ', 'integer',5};
defaults(end+1,:) = {'minHBratio', 'integer',.6};
defaults(end+1,:) = {'minAHratio', 'integer',.2};
defaults(end+1,:) = {'maxUncert', 'integer',2}; % pixels
defaults(end+1,:) = {'hybName', 'string',''}; % pixels
pars = ParseVariableArguments(varargin,defaults,mfilename); 

spotTable = CT4_BlankSpotTable(1);
xy2um = imProps.xy2um;
z2um = imProps.z2um;
        
zSteps = imProps.zSteps;
[~,fidMaxIdx] = max(cropIm(:));
[~,~,fidMax_z] = ind2sub(size(cropIm),fidMaxIdx);
zi = max(1,fidMax_z-pars.boxSizeZ);
ze = min(zSteps,fidMax_z+pars.boxSizeZ);

n=1; % loop over traces
t=1; % loop over traces

data_3d = cropIm(:,:,zi:ze);
% fit the 3D image
          fitTable = FitPsf3D(data_3d,...
                    'minPeakHeight',100,...
                    'maxFitWidth',pars.boxSize*2,...
                    'maxFitZdepth',pars.boxSize*2,...
                    'keepBrightest',1,...
                    'xyUnitConvert',1,...  108  % converted below
                    'zUnitConvert',1,... 250
                    'minHBratio',pars.minHBratio,...
                    'minAHratio',pars.minAHratio,...
                    'maxUncert',pars.maxUncert,...
                    'troubleshoot',false,...
                    'badSpots','nan');
            % convert some names, since we use the format written for the fit overlap spots below.  
            if isempty(fitTable) || isnan(fitTable.x)
                fitTable  = [];
                fitTable.x = nan;
                fitTable.y = nan;
                fitTable.z = nan;
                fitTable.z = nan;
                fitTable.b = nan;
                fitTable.a = nan;
                fitTable.wx = nan;
                fitTable.wz = nan;
                fitTable.fitSuccess = false;
            else
                fitTable.fitSuccess = true;
            end
            % convert some names
            fitTable.sXY = fitTable.wx;
            fitTable.sZ = fitTable.wz;
            fitTable.px = 0;
            fitTable.py = 0;                 

% "n" and "t" not needed, kept here for cross compatibility
spotTable.spotID(n) = nan; % will be populated in master table;
spotTable.traceID(n) = sptPixTable.traceID(t); %  
spotTable.x_um(n) = (-pars.boxSize + fitTable.x(1))*xy2um;
spotTable.y_um(n) = (-pars.boxSize + fitTable.y(1))*xy2um;
spotTable.z_um(n) = (-pars.boxSizeZ + fitTable.z(1))*z2um; 
spotTable.barcodeID(n) = nDataChns*(currHyb-1)+d;  % load from experiment table
spotTable.cellID(n) = sptPixTable.cellID(t);
% this is the additional data table  
spotTable.fitQuality(n) = fitTable.fitSuccess;
spotTable.x_fid(n) = sptPixTable.x_pix(t)*xy2um;
spotTable.y_fid(n) = sptPixTable.y_pix(t)*xy2um;
spotTable.z_fid(n) = (zi+pars.boxSizeZ)*z2um;  % center
%                 spotTable.x_shift_um(n) = shifts.xshift*xy2um;
%                 spotTable.y_shift_um(n) = shifts.yshift*xy2um;
%                 spotTable.z_shift_um(n) = shifts.zshift*z2um;
%                 spotTable.shift_score(n) = shifts.score;
 spotTable.spt_brightness(n) = fitTable.a(1);
%                 spotTable.x_um_ci95(n) = xErrWidth(nn)*xy2um;
%                 spotTable.y_um_ci95(n) = yErrWidth(nn)*xy2um;
%                 spotTable.z_um_ci95(n) = zErrWidth(nn)*z2um;
spotTable.spt_bkd(n) = fitTable.b;
spotTable.spt_bkd_xtilt(n) = fitTable.px;
spotTable.spt_bkd_ytilt(n) = fitTable.py;
spotTable.spt_sXY(n) = fitTable.sXY;
spotTable.spt_sZ(n) = fitTable.sZ;
spotTable.fid_brightness(n) = sptPixTable.brightness(t);
spotTable.dat_chn(n) = d;
% constant data
% these are actually constant for the same raw image-stack (same fov + hyb)  
% it is convienent to have these in the same table though. 
spotTable.stageX(n) = imProps.stageXY(1); % different stages may report different units (mm, um)? 
spotTable.stageY(n) = imProps.stageXY(2);
spotTable.fov(n) = sptPixTable.fov(t);%  f; 
spotTable.hyb(n) = currHyb;
% added
spotTable.hybName = pars.hybName;