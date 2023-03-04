function [fiducialAlignFrames,regData] = RegisterImagesFast(fiducialFrames,varargin)
% Registers a stack of images, provided as a cell array, to one of the
% images in the stack, denoted as the 'refHyb' (default is position 1). 
%
% Updates: 
% 2017-07-29, 
% Revised to return the translation rotation coordinates
%
% 2018-? 
% Revised into RegisterImages_light to accelerate for ChrTracer2
% 
% 2018-11-24 
% Revised into RegisterImagesFast to acclerate alignment for ChrTracer3 and
% to allow more robust rotation registration.
%
% Alistair Boettiger
% 2017-04-06
% 
% Copyright CC BY

global figureSavePath;

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
% key parameters
defaults(end+1,:) = {'previousAlignFrames', 'freeType',[]}; % can pass a list of previously aligned frames for speed
defaults(end+1,:) = {'alignContrastLow', 'fraction', .8}; % low image threshold for contrast balance prior to coarse alignment
defaults(end+1,:) = {'alignContrastHigh', 'fraction', .99}; % high threshold  for contrast balance prior to coarse alignment
defaults(end+1,:) = {'refHybe','integer',1}; % hybe to use to start alignment
defaults(end+1,:) = {'hybs','integer',inf}; % hybe to use to start alignment
% common fov parameters
defaults(end+1,:) = {'verbose', 'boolean', true}; 
defaults(end+1,:) = {'veryVerbose', 'boolean', false}; 
defaults(end+1,:) = {'showPlots', 'boolean', true}; 
defaults(end+1,:) = {'showExtraPlots', 'boolean', false}; 
% defaults(end+1,:) = {'saveData', 'boolean', false};  obsolete.
defaults(end+1,:) = {'saveFigs', 'boolean', true}; 
defaults(end+1,:) = {'stopOnError','boolean',false};
% CorrAlignFast Parameters
defaults(end+1,:) = {'maxSize', 'positive', 400}; % rescale all images to this size for alignment
defaults(end+1,:) = {'fineBox', 'freeType', []};  % perform fine scale alignment using a box of this size around the brightest point.
defaults(end+1,:) = {'fineUpsample', 'positive', 3};  
defaults(end+1,:) = {'maxShift', 'nonnegative', inf};
defaults(end+1,:) = {'gradMax', 'boolean', true};
defaults(end+1,:) = {'minGrad', 'float', -inf};
defaults(end+1,:) = {'angles','float',0}; % -10:1:10
defaults(end+1,:) = {'scales','float',1}; % -10:1:10
defaults(end+1,:) = {'fineMaxShift', 'nonnegative', 10};
defaults(end+1,:) = {'fineAngles','float',0}; % -1:.1:1
defaults(end+1,:) = {'fineScales','float',1}; % 0.95:0.01:.1.05
defaults(end+1,:) = {'fineCenter','array',[0,0]};
defaults(end+1,:) = {'verbose','boolean',false};
defaults(end+1,:) = {'showplot', 'boolean', true};
defaults(end+1,:) = {'fastDisplay', 'boolean', true};
defaults(end+1,:) = {'displayWidth', 'integer', 500};
defaults(end+1,:) = {'showExtraPlot', 'boolean', false};
defaults(end+1,:) = {'minFineImprovement', 'float', .5}; 
defaults(end+1,:) = {'showCorrAlign', 'boolean', false};
defaults(end+1,:) = {'fov', 'integer', []};

pars = ParseVariableArguments(varargin,defaults,mfilename);

if pars.saveFigs
   mkdir([figureSavePath,filesep,'CorrAlign'] );
end

numHybes = size(fiducialFrames,1);
if isempty(pars.previousAlignFrames)
    fiducialAlignFrames = cell(numHybes,1); 
else
    fiducialAlignFrames = pars.previousAlignFrames;
end


% initialize regData
% need to record nans so these convert to table properly
% I we could have made these a table from the beginning but it doesn't
% perserve backwards compatiability 
for h=1:numHybes
    regData(h).xshift = nan; %#ok<*AGROW>  % 
    regData(h).yshift = nan;
    regData(h).theta = nan;
    regData(h).rescale = nan;
    regData(h).xshift2 = nan;
    regData(h).yshift2 = nan;
    regData(h).theta2 = nan;
    regData(h).rescale2 = nan; 
end
h = pars.refHybe;
regData(h).xshift = 0;
regData(h).yshift = 0;
regData(h).theta = 0;
regData(h).rescale = 1;
regData(h).xshift2 = 0;
regData(h).yshift2 = 0;
regData(h).theta2 = 0;
regData(h).rescale2 = 1;

% loop over target hybes and align to previous or align to reference
h1 = max(fiducialFrames{pars.refHybe}(:,:,:),[],3);     % figure(1); clf; imagesc(h1);
h1 = IncreaseContrast(h1,'low',pars.alignContrastLow,'high',pars.alignContrastHigh); % the target frame

if pars.hybs == inf
    hybs = 1:numHybes;
else
    hybs = pars.hybs;
end


for h=hybs % loop over target hybes
    try
        if pars.verbose
            disp(['processing data from hyb ',num2str(h)]);
        end
       
        if  h==pars.refHybe
            fiducialAlignFrames{h}(:,:,:) = fiducialFrames{h}(:,:,:);
        elseif h~=pars.refHybe % h=15
            % Correct XY Drift using Z-max projected images   
            fidFrame_h = fiducialFrames{h}(:,:,:);
            h2 = max(fidFrame_h,[],3); % the current frame
            h2 = IncreaseContrast(h2,'low',pars.alignContrastLow,'high',pars.alignContrastHigh); 
            if pars.showCorrAlign
                f30 = figure(30); clf;
                showplot = true;
            else
                showplot = false;
            end
            pars.label1 = ['h',num2str(pars.refHybe)];
            pars.label2 = ['h',num2str(h)];
            if ~isempty(pars.fov)
                pars.label1 = ['fov',num2str(pars.fov),' ',pars.label1];
            end
            regData(h) = CorrAlignFast(h1,h2,'showplot',showplot,'parameters',pars); %#ok<AGROW>
            if pars.showCorrAlign
                if pars.saveFigs
                    SaveFigure(f30,'name',['CorrAlign',filesep,'CorrAlign_fov',num2str(pars.fov,'%04d'),'_h',num2str(h,'%04d')],'formats',{'png'},'verbose',pars.veryVerbose,'overwrite',true);
                end
                pause(.1); 
            end
            fiducialAlignFrames{h}(:,:,:) = ApplyReg(fidFrame_h,regData(h));
            if pars.veryVerbose
                disp(['rotating data by ',num2str(angle),' degrees']);
            end
        end
    catch er
        disp(er.getReport);
        warning(['error aligning data on hyb ',num2str(h)]);
    end

end

if pars.verbose
    disp('finished aligning data');
end