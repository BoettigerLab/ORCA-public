function sisRatio = DetectSisters(spotTrace_s,varargin)


defaults = cell(0,3); % 
% parameters used for pairing sisters between the color channels
defaults(end+1,:) = {'zDepth','integer',0}; % 0 for auto
pars = ParseVariableArguments(varargin,defaults,mfilename);

if pars.zDepth == 0
    zDepth = max(spotTrace_s(:,3,1));
else
    zDepth = pars.zDepth;
end
tObs = size(spotTrace_s,1);
tzObs = tObs/zDepth;

% re-organizing data into z-stacks for easier processing
    currStack1_s1 = nan(tzObs,8,zDepth); % each sister track
    currStack1_s2 = nan(tzObs,8,zDepth);
    for z =1:zDepth
        currStack1_s1(:,:,z) = squeeze(spotTrace_s(z:zDepth:end,:,1));
        currStack1_s2(:,:,z) = squeeze(spotTrace_s(z:zDepth:end,:,2));
    end
    d=1;
    obsPerZ_1 = squeeze(sum(~isnan(currStack1_s1(:,d,:)),1));
    obsPerZ_2 = squeeze(sum(~isnan(currStack1_s2(:,d,:)),1));
    tot_obs_s1 = sum(obsPerZ_1) ;
    tot_obs_s2 = sum(obsPerZ_2);
    zsplits = zeros(tzObs,zDepth-1);
    for z = 1:zDepth-1
        zsplits(:,z) =  sqrt(nansum( (currStack1_s1(:,1:2,z) - currStack1_s1(:,1:2,z+1)).^2,2)) > 2;
    end
    for z = 1:zDepth-2
        zsplits(:,z) = zsplits(:,z) +  (sqrt(nansum( (currStack1_s1(:,1:2,z) - currStack1_s1(:,1:2,z+2)).^2,2)) > 2);
    end
    tot_zsplits = sum(zsplits(:));
    tot_sister = tot_obs_s2 + tot_zsplits; 
    sisRatio = tot_sister / tot_obs_s1;
    % noSister = tot_sister < pars.minSisObs*tot_obs_s1;

    %%
% 
%  ts = 1:tObs;
% figure(3); clf; 
% ax = [];
% zsym = {'<','^','>','v','s'};
% for d=1:2
%     ax(d) = subplot(2,1,d);
%     for z=1:zDepth
%         plot(ts(z:zDepth:end),squeeze(spotTrace_s(z:zDepth:end,d,1)),['r',zsym{z}]); hold on;
%         plot(ts(z:zDepth:end),squeeze(spotTrace_s(z:zDepth:end,d,2)),['m',zsym{z}]);
%         md = nanmedian( squeeze(spotTrace_s(z:zDepth:end,d,1)) );
%         if z<zDepth
%             zsp = md.*zsplits(:,z); zsp(zsp==0)=nan;
%             plot(ts(z:zDepth:end),zsp ,['b',zsym{z}] );
%         end
%     end
% end
% linkaxes(ax,'x');%  title(['s=',num2str(s)])