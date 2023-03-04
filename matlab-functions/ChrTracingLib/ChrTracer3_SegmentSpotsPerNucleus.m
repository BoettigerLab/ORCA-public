function [pars] = ChrTracer3_SegmentSpotsPerNucleus(eTableXLS,varargin)
%  inputs -- just and eTableXLS
% CT3 loads data from disk, rather than passing from memory for each step.
% this facilitates resume from where we left off behavior
% 
% hint: set fov = inf and overwrite = true to auto reprocess all;
% 

%% parameters
defaults = cell(0,3);
defaults(end+1,:) = {'Tool',{'Cellpose','SpotSelector'},'Cellpose'}; % switch back to other tool
defaults(end+1,:) = {'dataType',{'fiducial','data','all'},'fiducial'}; % try using data if hybes don't work
defaults(end+1,:) = {'hybsToLoad', 'integer', 1};  % optionally, use multiple hybes. or All hybes from all data channels (if the fiducial is weak, for example)
defaults(end+1,:) = {'fov', 'integer', 1};  % 
defaults(end+1,:) = {'dataChannel', 'integer', inf}; % 
defaults(end+1,:) = {'laplace', 'boolean', true}; %  use laplace filter
defaults(end+1,:) = {'daxRootDefault', 'string', 'ConvZscan*.dax'}; 
defaults(end+1,:) = {'analysisFolder', 'string', ''}; % folder in which to find and fovNNN_regData.csv 
defaults(end+1,:) = {'fovFolder', 'string', ''}; % folder to save the cellpose input/ouput stuff. 
defaults(end+1,:) = {'simplifyOutput', 'boolean', false}; % convert cell-array to image matrix if only 1 movie is requested
defaults(end+1,:) = {'displayContrastHigh', 'fraction', .999}; % convert cell-array to image matrix if only 1 movie is requested
defaults(end+1,:) = {'displayContrastLow', 'fraction', .5}; % convert cell-array to image matrix if only 1 movie is requested
defaults(end+1,:) = {'rerunCellpose','boolean',false};
defaults(end+1,:) = {'overwrite','boolean',false};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'saveFig','boolean',true};
defaults(end+1,:) = {'nucleiThreshold','fraction',.985}; % used in cell segmentation, prior to passing to cellpose
defaults(end+1,:) = {'diameter','nonnegative',0}; % used in cell segmentation, 0 for auto
pars = ParseVariableArguments(varargin,defaults,mfilename);



% load data from Etable
daxData = LoadDaxFromEtable(eTableXLS,...
            'fov',pars.fov,...
            'hybNumber',pars.hybsToLoad,...
            'fixDrift',true,...
            'dataType',pars.dataType,...
            'daxRootDefault',pars.daxRootDefault,...
            'driftFolder',pars.analysisFolder,...
            'simplifyOutput',false,...
            'verbose',pars.verbose);

%% rewrite me so inf = autocycle
if isinf(pars.fov)
    fov = 1:size(daxData,2);
else
    fov = pars.fov;
end

for f=fov  
    imHybs = squeeze(daxData(:,f,:));
    nucIm = max(cat(3,imHybs{:}),[],3);  % some day we may want to allow a different image to be selected for cell segmentation

    % ----- select a particular channel for nuc image ---------
    % nucImage = LoadDaxFromEtable(eTableXLS,'fov',f,'hybNumber',10,'fixDrift',true,'dataType','fiducial','driftFolder',pars.analysisFolder,'maxProject',false); % select a hyb to use for segmentation
    % nucIm = nucImage(:,:,40);
    % figure(5); clf; imagesc(nucIm);
    % -------------------------------------------------------% 

    %----- Combine hybs and show spot image
    % (with or without laplace filter) 
    if pars.laplace
        nChns = length(imHybs);
        normIm = cell(nChns,1);
        for c=1:nChns
            temp = imHybs{c};
            temp = imfilter(double(temp),fspecial('laplacian'),'symmetric');  
            temp(temp>0) = 0; 
            temp = uint16(max(temp(:)) - temp); % figure(12); clf; imagesc(temp); colorbar;
            normIm{c} = temp;
        end
        % figure(1); clf; imagesc(Ncolor(25*cat(3,normIm{:})));  % just for troubleshooting 
        imSpots = mean(cat(3,normIm{:}),3);
    else
        if length(imHybs) > 1
            imSpots = max(cat(3,imHybs{:}),[],3);
        else
            imSpots = imHybs{1}; % this is faster than cat on nothing. 
        end
    end
    % imSpots = max(cat(3,normIm{:}),[],3);
    % imSpots(imSpots>quantile(imSpots(:),.999)) = max(imSpots(:));
    % figure(3); clf; imagesc(imSpots); caxis([20,2000])  
    dispSpot = IncreaseContrast(imSpots,'high',pars.displayContrastHigh,...
                                        'low',pars.displayContrastLow);
    dispNuc = IncreaseContrast(nucIm,'high',pars.displayContrastHigh,...
                                        'low',pars.displayContrastLow);
    figure(1); clf; Ncolor(cat(3,dispNuc,dispSpot)); colormap(gray); colorbar;  
        
      %% cellpose Segment 
    % if file exists, skip, unless overwrite is on
    fileName = [pars.analysisFolder,filesep,'fov',num2str(f,'%03d'),'_selectSpots.csv'];
    if exist(fileName,'file') && ~pars.overwrite
        % load and show
        spots = readtable(fileName);
        figure(1); hold on;
        plot(spots{:,1},spots{:,2},'yo');
        answer = questdlg('found existing fits, recompute?', ... % question
                'Prompt ', ...  % pop-up label
                'Yes','No','No'); % op1 op2 default
        if strcmp(answer,'Yes')
            do_fit = true;
        else
            do_fit = false;
            disp('advancing to next FOV');
            pars.fov = f+1;
        end
    else
        do_fit = true;
    end
    if do_fit
        fovFolder = [pars.analysisFolder,'fov',num2str(f,'%02d'),filesep];
        [fidTable,cellID_full,figH] = SegmentSpotsPerNucleus(imSpots,'f',f,...
                                        'nucImage',nucIm,...
                                        'saveFolder',fovFolder,...
                                        'overwrite',true,...
                                        'rerunCellpose',pars.rerunCellpose,...
                                        'nucleiThreshold',pars.nucleiThreshold);     
        %% approve results
        if ~isinf(pars.fov)
         answer = questdlg('save result or change parameters?', ... % question
                'Prompt ', ...  % pop-up label
                'Save and Continue','Refit','Save and Continue'); % op1 op2 default
        else
            answer ='Save and Continue';
        end

        if strcmp(answer,'Save and Continue')         
            %---- save results
            writetable(fidTable,fileName);
            if pars.verbose
                disp(['wrote ',fileName])
            end
            if pars.saveFig
                figure(figH);
                figName = [pars.analysisFolder,filesep,'fov',num2str(f,'%03d'),'_selectSpots'];
                % 2020a updated 
                savefig(figH,[figName,'.fig','compact']); 
                exportgraphics(figH,[figName,'.png'],'resolution',300); %  
                % Alternatively, save graphics,'fig.pdf','ContentType','vector')
            end
            disp('advancing to next FOV');
            if ~isinf(pars.fov)
                pars.fov = f+1;
            end
            % ---- update parameters and refit
            
         end



    end
end