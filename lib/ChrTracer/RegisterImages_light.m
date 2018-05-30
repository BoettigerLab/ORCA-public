function [fiducialAlignFrames,regData,goodHybes] = RegisterImages_light(fiducialFrames,varargin)
% 
%
% Updates: 
% 2017-07-29, 
% Revised to return the translation rotation coordinates
%
% Alistair Boettiger
% 2017-04-06
% 
% Copyright CC BY



% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
% key parameters
defaults(end+1,:) = {'alignmentBoxWidth', 'positive', inf}; % pixels.  size of region to use for frame alignment
defaults(end+1,:) = {'alignUpsample', 'positive', 1}; % upsample data for coarse alignment (slow, should be unnecessary);
defaults(end+1,:) = {'alignContrastLow', 'fraction', .8}; % low image threshold for contrast balance prior to coarse alignment
defaults(end+1,:) = {'alignContrastHigh', 'fraction', .99}; % high threshold  for contrast balance prior to coarse alignment
defaults(end+1,:) = {'refHybe','integer',1}; % hybe to use to start alignment
defaults(end+1,:) = {'previousAlignFrames','array',[]}; % previous align frames
defaults(end+1,:) = {'targetHybes','integer',[]};
defaults(end+1,:) = {'goodHybes','boolean',[]};
defaults(end+1,:) = {'corrAngles','freeType',[]}; % rotate angles to check for better alignment   
defaults(end+1,:) = {'rotation','boolean',true}; % also compute and correct rotation angle
defaults(end+1,:) = {'maxD','positive',15}; % distance in pixels to match objects prior to computing rotation
% common fov parameters
defaults(end+1,:) = {'verbose', 'boolean', true}; 
defaults(end+1,:) = {'veryVerbose', 'boolean', false}; 
defaults(end+1,:) = {'showPlots', 'boolean', true}; 
defaults(end+1,:) = {'showExtraPlots', 'boolean', false}; 
defaults(end+1,:) = {'saveData', 'boolean', false}; 
defaults(end+1,:) = {'stopOnError','boolean',false};

pars = ParseVariableArguments(varargin,defaults,mfilename);


numHybes = size(fiducialFrames,1);
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


% initialize regData
regData(1).xshift = 0;
regData(1).yshift = 0;
regData(1).angle = 0;
regData(1).tform = []; 
regData(numHybes).xshift = 0;
regData(numHybes).yshift = 0;
regData(numHybes).angle = 0;
regData(numHybes).tform = []; 

% loop over target hybes and align to previous or align to reference
h1 = max(fiducialFrames{pars.refHybe}(:,:,:),[],3);     % figure(1); clf; imagesc(h1);
h1 = IncreaseContrast(h1,'low',pars.alignContrastLow,'high',pars.alignContrastHigh); % the target frame


for h=pars.targetHybes % loop over target hybes
    try
        if pars.verbose
            disp(['processing data from hyb ',num2str(h)]);
        end
       
        if  h==pars.refHybe
            fiducialAlignFrames{h}(:,:,:) = fiducialFrames{h}(:,:,:);
        elseif h~=pars.refHybe % h=15
            % Correct Drift using max projected images             
            fidFrame_h = fiducialFrames{h}(:,:,:);
            h2 = max(fidFrame_h,[],3); % the current frame
            h2 = IncreaseContrast(h2,'low',pars.alignContrastLow,'high',pars.alignContrastHigh); 
            
            if pars.showExtraPlots;  figure(10); clf; end
            if ~isempty(pars.corrAngles)  % Correct arbitrary rotations (default is +/-10 degrees)
                [xshift,yshift,angle] = CorrAlignRotate(h1,h2,'region',pars.alignmentBoxWidth,'upsample',pars.alignUpsample,'showplot',pars.showExtraPlots,'angles',pars.corrAngles);
                if pars.veryVerbose
                    disp(['rotating data by ',num2str(angle),' degrees']);
                end
            else
                [xshift,yshift] = CorrAlign(h1,h2,'region',pars.alignmentBoxWidth,'upsample',pars.alignUpsample,'showplot',pars.showExtraPlots);
                % [xshift,yshift] = CorrAlign(h1(:,1:600),h2(:,1:600),'region',pars.alignmentBoxWidth,'upsample',pars.alignUpsample,'showplot',pars.showExtraPlots);
                % [xshift,yshift] = CorrAlign(h1(1:600,:),h2(1:600,:),'region',pars.alignmentBoxWidth,'upsample',pars.alignUpsample,'showplot',pars.showExtraPlots);
                % disp('aligning');
                angle = 0;
            end
            % Apply x-y shift and coarse rotation
            fidFrame_h = imrotate(ImageTranslate(fiducialFrames{h}(:,:,:),[xshift,yshift]),angle,'bilinear','crop');  % imtranslate
            
            
            % Correct Fine Scale Rotation by nearest neighbor mapping
            tform = []; 
            if pars.rotation 
                h2 = imrotate(ImageTranslate(h2,[xshift,yshift]),angle,'bilinear','crop');
                try
                    % tform = imregcorr(h2,h1,'rigid');         
                    % fidFrame_h = imwarp(fidFrame_h,tform,'OutputView',imref2d(size(fidFrame_h)));
                    
                    [~,tform] = RegisterRotatedImage(h1,h2,'maxD',pars.maxD,'showPlots',pars.showExtraPlots); % based on nearest neighbor matching of points
                    fidFrame_h = imwarp(fidFrame_h,tform.tobj,'OutputView',imref2d(size(fidFrame_h)));
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
            
            % save the shifts 
            regData(h).xshift = xshift;
            regData(h).yshift = yshift;
            regData(h).angle = angle;
            regData(h).tform = tform; 
            
            % save the aligned data
            fiducialAlignFrames{h}(:,:,:) = fidFrame_h;
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
