function [fiducialAlignFrames,dataAlignFrames,goodHybes,regData] = RegisterImages(fiducialFrames,varargin)
% 
%
% 
% Alistair Boettiger
% 2017-04-06
% Copyright CC BY

% 
if ~isempty(varargin)
    if ~ischar(varargin{1})
        dataFrames = varargin{1};
        varin = varargin(2:end);
    else
        dataFrames = fiducialFrames;
        varin = varargin;
    end
else
    dataFrames = fiducialFrames;
    varin = {};
end

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
% key parameters
defaults(end+1,:) = {'rotation','boolean',true}; % also compute and correct rotation angle
defaults(end+1,:) = {'maxD','positive',15}; % distance in pixels to match objects prior to computing rotation
defaults(end+1,:) = {'alignmentBoxWidth', 'positive', inf}; % pixels.  size of region to use for frame alignment
defaults(end+1,:) = {'alignUpsample', 'positive', 1}; % upsample data for coarse alignment (slow, should be unnecessary);
defaults(end+1,:) = {'alignContrastLow', 'fraction', .8}; % low image threshold for contrast balance prior to coarse alignment
defaults(end+1,:) = {'alignContrastHigh', 'fraction', .99}; % high threshold  for contrast balance prior to coarse alignment
defaults(end+1,:) = {'alignToFirst', 'boolean', true}; % align to first image? if false will align to previous non-empty image data
defaults(end+1,:) = {'minFracObj','fraction',.65}; % min fraction of objects from previous hybe found in target hybe to be acceptable 
defaults(end+1,:) = {'refHybe','integer',1}; % hybe to use to start alignment
defaults(end+1,:) = {'previousAlignFrames','array',[]}; % previous align frames
defaults(end+1,:) = {'targetHybes','integer',[]};
defaults(end+1,:) = {'goodHybes','boolean',[]};
defaults(end+1,:) = {'flattenBackground','nonnegative',0};
defaults(end+1,:) = {'corrAngles','freeType',[]}; % rotate angles to check for better alignment   
defaults(end+1,:) = {'subW','freeType',[]};
defaults(end+1,:) = {'subH','freeType',[]};
defaults(end+1,:) = {'warp',{'NonreflectiveSimilarity','lwm','polynomial','pwl'},'NonreflectiveSimilarity'};
% common fov parameters
defaults(end+1,:) = {'verbose', 'boolean', true}; 
defaults(end+1,:) = {'veryVerbose', 'boolean', false}; 
defaults(end+1,:) = {'showPlots', 'boolean', true}; 
defaults(end+1,:) = {'showExtraPlots', 'boolean', false}; 
defaults(end+1,:) = {'saveData', 'boolean', false}; 
defaults(end+1,:) = {'stopOnError','boolean',false};

pars = ParseVariableArguments(varin,defaults,mfilename);


[numHybes,numDataChns] = size(fiducialFrames);

if isempty(pars.targetHybes)
    pars.targetHybes = 1:numHybes;
end

if isempty(pars.goodHybes)
    pars.goodHybes = true(1,numHybes);
end

if isempty(pars.previousAlignFrames)
    fiducialAlignFrames = cell(numHybes,1); 
else
    fiducialAlignFrames = pars.previousAlignFrames;
end
if isempty(dataFrames)
    dataFrames = fiducialFrames;
end

dataAlignFrames = cell(numHybes,numDataChns);


% initialize count
% start a rough count of objects in image (used to determine image quality). 
h1 = max(fiducialFrames{pars.refHybe}(:,:,:),[],3);

% check if a sub-region was requested for correlation alignment
[h,w] = size(h1);
if isempty(pars.subH)
    subH = 1:h;
else
    subH = pars.subH;
end
if isempty(pars.subW)
    subW = 1:w;
else
    subW = pars.subW;
end


objMap = imregionalmax(imresize(h1,.02));
numObj = sum(objMap(:));

% initialize regData
regData(numHybes).xshift = 0;
regData(numHybes).yshift = 0;
regData(numHybes).angle = 0;
regData(numHybes).tform = []; 

% loop over target hybes and align to previous or align to reference
for h=pars.targetHybes % loop over target hybes
    try
        if pars.verbose
            disp(['processing data from hyb ',num2str(h)]);
        end
        if h~=pars.refHybe % h=15
            % Correct Drift using max projected images
            if pars.alignToFirst 
                alignHyb = pars.refHybe;
                h1 = max(fiducialAlignFrames{alignHyb}(:,:,:),[],3);
                if pars.flattenBackground > 0
                    bkd = imresize(imresize(h1,pars.flattenBackground),size(h1)); 
                    h1 = h1 - bkd;
                end
            else % alignToPrevious.  
                % If previous image is out of focus (substantial decrease
                % in the number of objects detected relative to image 1),
                % use the image before that one for alignment. 
                % This method is better for addressing the gradual change
                % that sometimes happens across the data set.  
                alignHyb = h-1; 
                numObjH = 0; 
                while numObjH < pars.minFracObj*numObj && alignHyb > 0
                    h1 = max(fiducialAlignFrames{alignHyb}(:,:,:),[],3);
                    objMap = imregionalmax(imresize(h1,.02));
                    numObjH = sum(objMap(:));
                    if numObjH < pars.minFracObj*numObj % record bad hybes; 
                        if pars.verbose
                           warning(['hybe ',num2str(alignHyb),' found only ',num2str(numObjH),' of ',num2str(numObj)]);
                        end
                        pars.goodHybes(alignHyb) = false; 
                    end
                    alignHyb = alignHyb -1;         
                end
            end

            h1 = IncreaseContrast(h1,'low',pars.alignContrastLow,'high',pars.alignContrastHigh); % the target frame
            h2 = max(fiducialFrames{h}(:,:,:),[],3); % the current frame
            if pars.flattenBackground > 0
                    bkd = imresize(imresize(h2,pars.flattenBackground),size(h2)); 
                    h2 = h2 - bkd;
            end
            h2 = IncreaseContrast(h2,'low',pars.alignContrastLow,'high',pars.alignContrastHigh); 
            if pars.showExtraPlots; corrFig = figure(10); clf; end
            if ~isempty(pars.corrAngles)  % Correct arbitrary rotations (default is +/-10 degrees)
                [xshift,yshift,angle] = CorrAlignRotate(h1(subH,subW),h2(subH,subW),'region',pars.alignmentBoxWidth,'upsample',pars.alignUpsample,'showplot',pars.showExtraPlots,'angles',pars.corrAngles);
                if pars.veryVerbose
                    disp(['rotating data by ',num2str(angle),' degrees']);
                end
            else
                [xshift,yshift] = CorrAlign(h1(subH,subW),h2(subH,subW),'region',pars.alignmentBoxWidth,'upsample',pars.alignUpsample,'showplot',pars.showExtraPlots);
                angle = 0;
            end
            % Apply x-y shift and coarse rotation
            fidFrame_h = imrotate(ImageTranslate(fiducialFrames{h}(:,:,:),[xshift,yshift]),angle,'bilinear','crop');  % imtranslate
            datFrame_h = cell(1,numDataChns);
            for n=1:numDataChns
                datFrame_h{n} = imrotate( ImageTranslate(dataFrames{h,n}(:,:,:),[xshift,yshift]), angle,'bilinear','crop');
            end
            tform = []; 
            if pars.rotation % Correct Fine Scale Rotation
                h2 = imrotate(ImageTranslate(h2,[xshift,yshift]),angle,'bilinear','crop');
                try
                    [~,tform] = RegisterRotatedImage(h1,h2,'maxD',pars.maxD,'showPlots',pars.showExtraPlots,'warp',pars.warp); % based on nearest neighbor matching of points
                    fidFrame_h = imwarp(fidFrame_h,tform.tobj,'OutputView',imref2d(size(fidFrame_h)));
                    for n=1:numDataChns
                        datFrame_h{n} = imwarp(datFrame_h{n},tform.tobj,'OutputView',imref2d(size(datFrame_h{n})));
                    end
                    if pars.veryVerbose
                        disp(['rotated data by ',num2str(tform.rotationAngle),' degrees']);
                    end
                catch er
                    if pars.veryVerbose
                        disp(er.message);
                        disp(['failed to rotate hybe ',num2str(h)]);
                    end
                end              
            end
            regData(h).xshift = xshift;
            regData(h).yshift = yshift;
            regData(h).angle = angle;
            regData(h).tform = tform; 
            fiducialAlignFrames{h}(:,:,:) = fidFrame_h;
            for n=1:numDataChns
                dataAlignFrames{h,n}(:,:,:) = datFrame_h{n};  
            end
        else
            h1 = IncreaseContrast(h1,'low',pars.alignContrastLow,'high',pars.alignContrastHigh); % the target frame
            fiducialAlignFrames{h}(:,:,:) = fiducialFrames{h}(:,:,:);
            for n=1:numDataChns
                dataAlignFrames{h}(:,:,:) = dataFrames{h}(:,:,:);
                dataAlignFrames{h,n}(:,:,:) = dataFrames{h,n}(:,:,:);
            end
        end
    catch er
        if pars.verbose; disp(er.getReport); end
        warning(['error aligning data on hyb ',num2str(h)]);
    end

end

if pars.verbose
    disp('finished aligning data');
end

goodHybes = pars.goodHybes;
