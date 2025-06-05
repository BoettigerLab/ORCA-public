function [currTrace1,currTrace2,unlinked1,unlinked2] = LinkZpairs(currStack1,currStack2,varargin)
%
% Approach: link traces across the z-stack if the frame-to-frame step-size
% difference per point after the merge is on the same scale as before the
% merge.  Same scale is defined by the maxFoldChange (2-fold default) and
% an absolute max size (1000 nm). 
%  The z-position with the most *shared* observations between the two
%  channels is used as the seed position. This reduces the chance of
%  picking different sisters for the two channels. 
%         Other positions in the stack are merged into this or
% rejected based on whether they are consistent with the path or behave
% like independent particles. 
% Second loci (such as on the sister chromosome) are excluded based on the 
% observation that they follow a different path from the original seed. If
% they follow the same path, they are intepreted as the same original
% allele seen from consecutive z-planes. 
% 
% 
% 
%% Inputs  
% currStack1
% currStack2  - two 3D matrix, t_Obs x d_Dims x z_Depth.   
%             d-Dims are x,y,z,t,h,b: x-position, y-position, z-position,
%             time/frame number, gaussian-fit height and gaussian fit
%             background.  
%             t-Obs is number of total z-stacks (not total frames)
%% Outputs
% currTrace1
% currTrace2 - a 2D matrix t_Time x d_Dims.  t_Time is now per frame
%            t_Time = t_Obs x z_Depth.
% 
%% updates
%  evolved from LinkZ.m, this version uses both traces 


%% optional pars
defaults = cell(0,3);
defaults(end+1,:) = {'maxFoldChange','nonnegative',2}; % distance in fold change to accept mergeing a localization 
defaults(end+1,:) = {'maxStepSize','nonnegative',1}; % distance in microns to accept merging a localization 
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'showFig','integer',10}; % note, uses 3 figures, 
defaults(end+1,:) = {'spotNum','integer',nan};
defaults(end+1,:) = {'npp','positive',108}; % nm per pixel (just for plotting/ troubleshooting)
pars = ParseVariableArguments(varargin,defaults,mfilename);

%% more parameters 
% 
% % data for troubleshooting
% currStack = spotStackCell{8};
% pars = ParseVariableArguments([],defaults,mfilename);

% a few more parameters translated from the inputs
[tObs,dDims,zDepth] = size(currStack1); 
s= pars.spotNum; 


%% main function

% find the longest z-sequence
ObsPerZ_1 = squeeze(sum(sum(~isnan(currStack1),2),1));
ObsPerZ_2 = squeeze(sum(sum(~isnan(currStack2),2),1));
ObsPerZ_Shared = squeeze(sum(sum(  (~isnan(currStack1) & ~isnan(currStack2)) ,2),1));

ObsPerZ = [ObsPerZ_1,ObsPerZ_2,ObsPerZ_Shared];

[~,zOrder] = sort(ObsPerZ_Shared,'descend');

currStacks = {currStack1,currStack2};
for chn=1:2
    sortStack = currStacks{chn}(:,:,zOrder);
    % sortStack(200,:,1)
    % sortStack(200,:,2)
    currTrace = nan(tObs*zDepth,dDims);
    unlinked = nan(tObs*zDepth,dDims);
    % we start using the data from the z-plane with the most data
    %   (after sorting, this is row 1)
    currTrace(zOrder(1):zDepth:end,:) = sortStack(:,:,1); 
    %  we compute the variation of the longest unmerged trace. This is our
    %  baseline stepsize for the particle.  If we find the same particle in
    %  other frames, it should fit within this trajectory by this amount. 
    stpVar = nanmean(sqrt( (diff(sortStack(:,1,1))).^2 + (diff(sortStack(:,2,1))).^2   ));
    % for future: might make more sense to use a "local" step size instead of a
    % global one.  or even just the previous step. 
    
    
    if pars.showFig
        figure(pars.showFig); clf;
    end
    for z=1:zDepth-1 % z=4
        % now we combine the next z-plane with this one;  
        traceMerge =  nan(tObs*2,dDims); % an array to allow us to compare the step size differences using the diff function 
        traceMerge(1:2:end,:) =sortStack(:,:,1); % this is our base data
        traceMerge(2:2:end,:) = sortStack(:,:,z+1); % each additional z in turn is considered as an addition to the first z
        mergeVar = sqrt( (diff(traceMerge(:,1))).^2 + (diff(traceMerge(:,2))).^2   )';
        accept_stepUp = [false,mergeVar < pars.maxFoldChange*stpVar  | mergeVar < pars.maxStepSize];
        accept_stepDown = [mergeVar < pars.maxFoldChange*stpVar  | mergeVar < pars.maxStepSize,false];
        accept = accept_stepUp | accept_stepDown;
        accept(1:2:end) = true; % always keep data from starting / longest trace 
        traceMerge0 = traceMerge;
        notLinked = traceMerge;
        notLinked(1:2:end,:) = nan; % reference plane is always linked 
        % notLinked(accept,:) = nan;  % data added from other planes are removed (might skip this, its okay to share)   
        traceMerge(~accept,:) = nan;
        currTrace(zOrder(z+1):zDepth:end,:) = traceMerge(2:2:end,:);
        unlinked(zOrder(z+1):zDepth:end,:) = notLinked(2:2:end,:);
           
        if pars.showFig
            % --- Just plotting
            x1 = z:zDepth:(tObs*zDepth);
            x2 = (z+1):zDepth:(tObs*zDepth);
            xMerge = nan(1,tObs*2);
            xMerge(1:2:end) = x1;
            xMerge(2:2:end) = x2;
            xx = xMerge;
            figure(pars.showFig+1); clf; 
            axx1 = subplot(2,1,1);
            plot(x1,sortStack(:,1,1),'.-','linewidth',1); hold on;
            plot(x2,sortStack(:,1,z+1),'.-'); hold on;
            axx2 = subplot(2,1,2);
            plot(xMerge,traceMerge0(:,1),'.-','linewidth',1); hold on;
            yy  = traceMerge(:,1);
            plot(xx(~isnan(yy)),yy(~isnan(yy)),'.-'); hold on;
            linkaxes([axx1,axx2],'x');
            legend();
        
            figure(pars.showFig); 
            x = z:zDepth:tObs*zDepth;
            ax1 = subplot(3,1,1); plot(x,squeeze(sortStack(:,1,z)),'.-'); hold on; title(s);
            ax2 = subplot(3,1,2); plot(x,squeeze(sortStack(:,2,z)),'.-'); hold on; title(s);
            ax3 = subplot(3,1,3); plot(x,squeeze(sortStack(:,5,z)),'.-'); hold on; title(s);
            linkaxes([ax1,ax2,ax3],'x');
            % does the distance between the points balance out the z-associated
            % step bias.  
        end
    end
    
    
    if pars.showFig
        figure(pars.showFig); legend();
        xx = 1:(tObs*zDepth);
        yy = currTrace(:,1);  % just the x-coord for simplicity   
        figure(pars.showFig + 2); clf; 
        plot(xx(~isnan(yy)),yy(~isnan(yy)),'.-'); hold on;
        ys = sortStack(:,1,1);
        orig_aveXstep = nanmedian(abs(diff(ys))) ;
        aveXstep = nanmedian(abs(diff(yy))) ;
        fracGained = sum(~isnan(yy))./sum(~isnan(ys));
        fracSkipped = 1-sum(ObsPerZ)/max(ObsPerZ);
    
         title({['Spot ',num2str(pars.spotNum)]; [' x-step per frame: '  ,num2str(pars.npp*aveXstep,3), ' nm. '...
            '  orig per ', num2str(zDepth), ' frames: ',num2str(pars.npp*orig_aveXstep,3), ' nm. '];...
            [' Skipped fraction=',num2str(fracSkipped,3) ,'.   Coverage increase= ' num2str(fracGained,3)]});   
        % note - if no consecutive spots are found, the average step-to-step
        % variation will be NaN.  
    end
    % export
    if chn==1
        currTrace1 = currTrace;
        unlinked1 = unlinked;
    elseif chn==2
        currTrace2 = currTrace;
        unlinked2 = unlinked;
    end
end

