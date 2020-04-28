function [fidSpots,dataSpots,pars] = ChrTracer2_CropAllSpots(fidMapData,datMapData,varargin)
% pars = ChrTracer_CropAndPlot(fidMapData,datMapData,'parameters',pars);
% pars = ChrTracer_CropAndPlot(fidMapData,datMapData,'parameters',pars,'optionName',optionName);
% 
% -------------------------------------------------------------------------
% Required Inputs
% -------------------------------------------------------------------------
% Memory map files for fiducial and data: 'fidMapData', 'datMapData'
% 'parameters',pars
% 
% -------------------------------------------------------------------------
% Optional inputs
% -------------------------------------------------------------------------
% 
% 
% -------------------------------------------------------------------------
% Outputs
% -------------------------------------------------------------------------
% 
% 
% -------------------------------------------------------------------------
% Notes
% -------------------------------------------------------------------------
% called by ChrTracer

% supress some unnecessary warnings. 
warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');
warning('off','MATLAB:prnRenderer:opengl');

defaults = cell(0,3);
% key parameters
defaults(end+1,:) = {'lociXY', 'array', []}; % actually a REQUIRED variable 
defaults(end+1,:) = {'boxWidth', 'positive', 16};
defaults(end+1,:) = {'goodHybes','array',[]};
defaults(end+1,:) = {'showFolderNames', 'boolean',false};
% FOV parameters
defaults(end+1,:) = {'fov', 'integer', 0};  % Field of view in experiment,  0 for not known;
defaults(end+1,:) = {'showPlots','boolean',true};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'selectSpots','integer',[]};
defaults(end+1,:) = {'numParallel','integer',1};
pars = ParseVariableArguments(varargin, defaults, mfilename);



%%
disp('Cropping Spots...');
% A little more parameter parsing

% Get some data sizes
[numHybes,numDataChns] = size(datMapData); 
numSpots = size(pars.lociXY,1);

% parse select spots
if isempty(pars.selectSpots)
    selectSpots = 1:numSpots;
else
    selectSpots = pars.selectSpots;
end

% determine if any data were flagged as to skip
pars.goodHybes = logical(pars.goodHybes);
if isempty(pars.goodHybes)
    pars.goodHybes = true(1,numHybes);
end


%% Crop images
tic;
try
    boxWidth = pars.boxWidth; % image width
    spots = pars.lociXY;
    fidSpots = cell(numSpots,1);
    dataSpots = cell(numSpots,1);
    im = ReadFromMemMap(fidMapData{1});
    [numRows,numCols,~] = size(im);


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

    if pars.numParallel == 1
        for s=selectSpots    
            yi = max(1,spots(s,2)-boxWidth/2);
            ye = min(spots(s,2)+boxWidth/2,numRows);
            xi =  max(1,spots(s,1)-boxWidth/2);
            xe = min(spots(s,1)+boxWidth/2,numCols);
            for h=find(pars.goodHybes)
                % load drift and rotation corrected image;
                % rotation must be applied to whole image, not to part. 
                fidSpots{s}{h} =  ReadFromMemMap(fidMapData{h},'roi',[xi,xe,yi,ye]); %  figure(10); clf; imagesc(max(temp,[],3));
                for n=1:numDataChns
                    dataSpots{s}{h,n} =ReadFromMemMap(datMapData{h,n},'roi',[xi,xe,yi,ye]); 
                end
            end
        end
    else   
        goodHybes = pars.goodHybes;
        parfor s=selectSpots    
            yi = max(1,spots(s,2)-boxWidth/2);
            ye = min(spots(s,2)+boxWidth/2,numRows);
            xi =  max(1,spots(s,1)-boxWidth/2);
            xe = min(spots(s,1)+boxWidth/2,numCols);
            for h=find(goodHybes)
                % load drift and rotation corrected image;
                % rotation must be applied to whole image, not to part. 
                fidSpots{s}{h} =  ReadFromMemMap(fidMapData{h},'roi',[xi,xe,yi,ye]); %  figure(10); clf; imagesc(max(temp,[],3));
                for n=1:numDataChns
                    dataSpots{s}{h,n} =ReadFromMemMap(datMapData{h,n},'roi',[xi,xe,yi,ye]); 
                end
            end
        end
    end

catch er
    disp(['encountered an error after  t=',num2str(t/60),'min.']);
    warning(er.getReport);
    error(er.message);
end

t = toc;
if pars.verbose
    disp(['finished crop and plot in t=',num2str(t/60),'min.']);
end
