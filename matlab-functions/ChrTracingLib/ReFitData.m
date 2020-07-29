function [dataTablesOut,allNewShifts] = ReFitData(fidS,datS,varargin)
%% recompute fine scale fiducial alignment and spot-fitting
% 
% Just for a specific spot in field of view. 
% datS = datSpots(:,:,:,:,s); 
% fidS = fidSpots(:,:,:,:,s);

defaults = cell(0,3);
% drift correction
defaults(end+1,:) = {'upsample','positive',4};
defaults(end+1,:) = {'maxXYdrift','positive',10};
defaults(end+1,:) = {'maxZdrift','positive',10};
defaults(end+1,:) = {'showExtraPlots','boolean',false};
defaults(end+1,:) = {'veryverbose','boolean',false};
% fitting 
defaults(end+1,:) = {'nmXYpix','positive',nmXY};
defaults(end+1,:) = {'nmZpix','positive',nmZ};
defaults(end+1,:) = {'maxXYstep','positive',12}; % max distance from fiducial center (first trim)
defaults(end+1,:) = {'maxZstep','positive',12}; % max distance from fiducial center (first trim)
defaults(end+1,:) = {'datMinPeakHeight', 'positive', 0};
defaults(end+1,:) = {'datTroubleshoot', 'boolean', false};
defaults(end+1,:) = {'datMaxFitWidth', 'positive', 12}; % Total number of pixels used for fitting, box centered on the brightest pixel in the trimmed area (2x maxXYstep for no further trimming)
defaults(end+1,:) = {'datMaxFitZdepth', 'positive', 12};
defaults(end+1,:) = {'datMinHBratio','nonnegative',.5}; %1.2 peak value over background value
defaults(end+1,:) = {'datMinAHratio','nonnegative',.1}; %0.25 fitted height over background vs peak value
defaults(end+1,:) = {'datMaxUncert','nonnegative',3}; %2 pixels
defaults(end+1,:) = {'bkdFrac','fraction',.1}; % 0 for off. dimmest bkdFrac of images will be used to compute local illumination background 
pars = ParseVariableArguments([],defaults,'fitPsf3D');
    
% =============== Compute 3D fine-scale drift correction ===============
tic
newShifts = cell(nH,1);
refHyb = 1;
refFid = fidS(:,:,:,refHyb);
fidTable = FindPeaks3D(refFid,...
                        'maxFitWidth',12,...
                        'keepBrightest',1,...
                        'troubleshoot',false);
for h=1:nH
[~,newShifts{h}] = Register3D(refFid,fidS(:,:,:,h),...
                'center',[fidTable.x(1),fidTable.y(1),fidTable.z(1)],...
                'upsample',pars.upsample,'verbose',pars.veryverbose,...
                'maxShiftXY',pars.maxXYdrift,'maxShiftZ',pars.maxZdrift,...
                'showplots',pars.showExtraPlots,...
                'returnImage',false); % just for debugging
end
allNewShifts = cellfun(@struct2table,newShifts,'UniformOutput',false);
allNewShifts = cat(1,allNewShifts{:});
toc
%=====================================================================%

% ================= Re-do data fitting ===============================
tic
% Data fitting defaults
fidXYZ = [fidTable.x(1),fidTable.y(1),fidTable.z(1)];
dataTables = cell(nH,1);

selHybes = 1:nH; % nH
for h=selHybes   
    xi = max(round(fidXYZ(1)-pars.maxXYstep ) ,1)  ;  % +allNewShifts.xshift(h)  (doesnt do anything inside the max)
    xe = min(round(fidXYZ(1)+pars.maxXYstep ) ,nX);
    yi = max(round(fidXYZ(2)-pars.maxXYstep ) ,1);
    ye = min(round(fidXYZ(2)+pars.maxXYstep ) ,nY);
    zi = max(round(fidXYZ(3)-pars.maxZstep  ) ,1);
    ze = min(round(fidXYZ(3)+pars.maxZstep  ) ,nZ);
    datSpt = datS(yi:ye,xi:xe,zi:ze,h); 
%     figure(5); clf; imagesc( max(datSpt,[],3));
%     figure(10); clf;
    dTable = FitPsf3D(datSpt,...
                    'minPeakHeight',pars.datMinPeakHeight,...
                    'maxFitWidth',pars.datMaxFitWidth,...
                    'maxFitZdepth',pars.datMaxFitZdepth,...
                    'keepBrightest',1,...
                    'xyUnitConvert',pars.nmXYpix,...
                    'zUnitConvert',pars.nmZpix,...
                    'minHBratio',pars.datMinHBratio,...
                    'minAHratio',pars.datMinAHratio,...
                    'maxUncert',pars.datMaxUncert,...
                    'troubleshoot',pars.datTroubleshoot);  

    % positions (add back the trim)
    if ~isempty(dTable)
        dTable.x = dTable.x + (xi-1)*pars.nmXYpix    +allNewShifts.xshift(h)*pars.nmXYpix ;  % in nm
        dTable.y = dTable.y + (yi-1)*pars.nmXYpix    +allNewShifts.yshift(h)*pars.nmXYpix;  % in nm
        dTable.z = dTable.z + (zi-1)*pars.nmZpix     +allNewShifts.zshift(h)*pars.nmZpix; % in nm
        dTable.xL = dTable.xL + (xi-1)*pars.nmXYpix  +allNewShifts.xshift(h)*pars.nmXYpix;  % in nm
        dTable.yL = dTable.yL + (yi-1)*pars.nmXYpix  +allNewShifts.yshift(h)*pars.nmXYpix; % in nm
        dTable.zL = dTable.zL + (zi-1)*pars.nmZpix   +allNewShifts.zshift(h)*pars.nmZpix;% in nm
        dTable.xU = dTable.xU + (xi-1)*pars.nmXYpix  +allNewShifts.xshift(h)*pars.nmXYpix;  % in nm
        dTable.yU = dTable.yU + (yi-1)*pars.nmXYpix  +allNewShifts.yshift(h)*pars.nmXYpix;  % in nm
        dTable.zU = dTable.zU + (zi-1)*pars.nmZpix   +allNewShifts.zshift(h)*pars.nmZpix;% in nm
        dTable.xshift = allNewShifts.xshift(h)*pars.nmXYpix;
        dTable.yshift = allNewShifts.yshift(h)*pars.nmXYpix;
        dTable.zshift = allNewShifts.zshift(h)*pars.nmZpix;
        dTable.hybe = h;
        dataTables{h} = dTable;
        
%         figure(5); hold on; 
%         plot(dTable.x,dTable.y,'r+');
%         title(num2str(h));
       %  pause;
    end
end  
toc
dataTablesOut = cat(1,dataTables{:});


