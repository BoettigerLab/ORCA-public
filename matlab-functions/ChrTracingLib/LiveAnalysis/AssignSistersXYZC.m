function [spotTrace1_out_s,spotTrace2_out_s] = AssignSistersXYZC(spotTrace1_s,spotTrace2_s,varargin) 
% Uses data from all sources: xy-sister doublets, z-stack doublets, and
% best correlation between data channels, simultaneously, in order to find
% the most parsamonious sister assignments. 
% 
% 
% AssignSisters runs on individual sisters, and doesn't guarentee that the
% which sister in chn1 matches which sister in channel 2. Moreover, when
% the sisters collide, it risks swapping their identies which is
% problematic of the alterant channel doesn't swap identities as well to
% keep the pair.  
% After calling AssignSister, this function considers the collisions, and
% consults the cognate color channel to decide if there has been an
% erroneous swap or not. 

% this function also doesn't need the whole spotTrace stack...


defaults = cell(0,3); % 
% parameters used for pairing sisters between the color channels
defaults(end+1,:) = {'zDepth','integer',0}; % 0 for auto
defaults(end+1,:) = {'verbose','boolean',false};
defaults(end+1,:) = {'swapDist','nonnegative',3}; % % entertains the potential to swap sister references paths that come within this distance (in xy-pixels) of one another. 
defaults(end+1,:) = {'minSister','nonnegative',1.5}; % 
defaults(end+1,:) = {'maxDistToRef','nonnegative',5}; % 
% parameters passed to AssignSisters
defaults(end+1,:) = {'minSisObs','nonnegative',0}; % 
defaults(end+1,:) = {'minSisWindows','nonnegative',6}; % 
defaults(end+1,:) = {'movAveWin','positive',50}; % 
defaults(end+1,:) = {'minObsPerAve','positive',4}; % min number of observations in moving average window to count (per chn, per z)
defaults(end+1,:) = {'subSample','positive',2}; % sampling interval per moving average window (values must be from 1 to size moveAveWin)
defaults(end+1,:) = {'maxDistToSeed','positive',30}; % 
defaults(end+1,:) = {'maxDistToPrev','positive',8}; %
defaults(end+1,:) = {'maxSisterDist','positive',45}; %
defaults(end+1,:) = {'figSis','nonnegative',0};
defaults(end+1,:) = {'figMovAve','nonnegative',0};
defaults(end+1,:) = {'figSisLabel','nonnegative',5};
defaults(end+1,:) = {'saveFig','boolean',true};
defaults(end+1,:) = {'traceNum','integer',0};
pars = ParseVariableArguments(varargin,defaults,mfilename);


swapDist = pars.swapDist; % entertains the potential to swap sister references paths that come within this distance (in xy-pixels) of one another. 
            % note, sister reference paths are downsampled to be smooth, so
            % they don't hold to 1-pixel of the data. 
minSister = pars.minSister; % sister reference lines closer than this will be considered singlets
maxDistToRef = pars.maxDistToRef; % points further than this distance in pixels will be excluded 

 %  92;% 33;  16  ~102

% 1 - singlet, loses all the data
% 52 - loses all the data
% 2 - tosses most of the data to the sister channel. 
s= pars.traceNum; 

% s=33; 31, 92
[chn1_ref_sis,chn1_conf] = AssignSisters( spotTrace1_s ,...
    'movAveWin',pars.movAveWin,'minObsPerAve',pars.minObsPerAve,'maxDistToSeed',pars.maxDistToSeed,'maxDistToPrev',pars.maxDistToPrev,'maxSisterDist',pars.maxSisterDist,...
    'minSisObs',pars.minSisObs,'figSis',pars.figSis);

sisChns1 = nansum(chn1_ref_sis(:,:,1),2) > 0;
if sisChns1(2) == 0
    chn1_ref_sis(2,:,:) = inf;
end

[chn2_ref_sis,chn2_conf] = AssignSisters(squeeze(  spotTrace2_s ),...
    'movAveWin',pars.movAveWin,'minObsPerAve',pars.minObsPerAve,'maxDistToSeed',pars.maxDistToSeed,'maxDistToPrev',pars.maxDistToPrev,'maxSisterDist',pars.maxSisterDist,...
    'minSisObs',0,'figSis',pars.figSis);

sisChns2 = nansum(chn2_ref_sis(:,:,1),2) > 0;
if sisChns2(2) == 0
    chn2_ref_sis(2,:,:) = inf;
end

if pars.zDepth == 0
    zDepth = max(spotTrace1_s(:,3,1));
else
    zDepth = pars.zDepth;
end
tObs = size(spotTrace1_s,1);
tzObs = tObs/zDepth;

% assume if you start closer to one sister you finish closer to that sister
d11 =  squeeze(sqrt(nansum( (chn1_ref_sis(1,:,:) - chn2_ref_sis(1,:,:)).^2 ,3)));
d12 =  squeeze(sqrt(nansum( (chn1_ref_sis(1,:,:) - chn2_ref_sis(2,:,:)).^2 ,3)));
d21 =  squeeze(sqrt(nansum( (chn1_ref_sis(2,:,:) - chn2_ref_sis(1,:,:)).^2 ,3)));
d22 =  squeeze(sqrt(nansum( (chn1_ref_sis(2,:,:) - chn2_ref_sis(2,:,:)).^2 ,3)));
dchn1 = squeeze(sqrt(nansum( (chn1_ref_sis(1,:,:) - chn1_ref_sis(2,:,:)).^2 ,3))) ;
dchn2 = squeeze(sqrt(nansum( (chn2_ref_sis(1,:,:) - chn2_ref_sis(2,:,:)).^2 ,3))) ;
% 
% figure(3); clf;
% subplot(4,1,1); 
%     d=1;
%     plot(chn1_ref_sis(1,:,d),'b.-'); hold on;
%     plot(chn1_ref_sis(2,:,d),'c.-');
%     plot(chn2_ref_sis(1,:,d),'r.-');
%     plot(chn2_ref_sis(2,:,d),'m.-');
%     title('x trace')
%     subplot(4,1,2); 
%     d=2;
%     plot(chn1_ref_sis(1,:,d),'b.-'); hold on;
%     plot(chn1_ref_sis(2,:,d),'c.-');
%     plot(chn2_ref_sis(1,:,d),'r.-');
%     plot(chn2_ref_sis(2,:,d),'m.-');
%     title('y trace')
% subplot(4,1,3);
%     title('sister distance')
%     plot(d11,'.-'); hold on;
%     plot(d12,'.-')
%     plot(dchn1< 1,'m+');
%      plot(dchn2< 1,'cx'); 
% subplot(4,1,4);
%     title('sister distance')
%     plot(d21,'.-'); hold on;
%     plot(d22,'.-')
%     plot(dchn1< 1,'m+');
%     plot(dchn2< 1,'cx'); 


% swap so ch1 s1 is paired with ch2 s1 if not already
kep = d11 + d22 < d12 + d21; % d11 < d12 & d22 < d21;
swp = d12 + d21 < d11 + d22; %  d12 < d11 & d21 < d22;  % 
ref1 = chn1_ref_sis; % short hand
ref2 = chn2_ref_sis;
ref1_conf = chn1_conf;
ref2_conf = chn2_conf;
if sum(swp) > sum(kep)
    ref1(1,:,:) = chn1_ref_sis(2,:,:);
    ref1(2,:,:) = chn1_ref_sis(1,:,:);
    ref1_conf(1,:,:) = chn1_conf(2,:,:);
    ref1_conf(2,:,:) = chn1_conf(1,:,:);
    if pars.verbose
        disp('swapping channels')
    end
    % only swap 1, otherwise this is just circular
end
% ref1a = ref1;
% ref2a = ref2;

% % recompute distances (for sanity in labels)
% d11 =  squeeze(sqrt(nansum( (ref1(1,:,:) - ref2(1,:,:)).^2 ,3)));
% d12 =  squeeze(sqrt(nansum( (ref1(1,:,:) - ref2(2,:,:)).^2 ,3)));
% d21 =  squeeze(sqrt(nansum( (ref1(2,:,:) - ref2(1,:,:)).^2 ,3)));
% d22 =  squeeze(sqrt(nansum( (ref1(2,:,:) - ref2(2,:,:)).^2 ,3)));
% dchn1 = squeeze(sqrt(nansum( (ref1(1,:,:) - ref1(2,:,:)).^2 ,3))) ;
% dchn2 = squeeze(sqrt(nansum( (ref2(1,:,:) - ref2(2,:,:)).^2 ,3))) ;

% find all the places where the sisters collide
posSwaps1 = regionprops(dchn1 < swapDist,'Centroid');
posSwaps1 = cat(1,posSwaps1.Centroid);
if ~isempty(posSwaps1)
    posSwaps1 = round(posSwaps1(:,1));
end
posSwaps2 = regionprops(dchn2 < swapDist,'Centroid');
posSwaps2 = cat(1,posSwaps2.Centroid);
if ~isempty(posSwaps2)
    posSwaps2 = round(posSwaps2(:,1));
end

if ~isempty(posSwaps1)
    % if swapping sister identity at the trace
    posSwaps1 = [0; posSwaps1; tObs];
    for p=1:length(posSwaps1)-1
        tw = posSwaps1(p)+1:posSwaps1(p+1);
        % compute distance
        r11 = ref1(1,tw,:);
        r21 = ref2(1,tw,:);
        r12 = ref1(2,tw,:);
        r22 = ref2(2,tw,:);
        ds11 =  nanmean( squeeze(sqrt(nansum( (r11 - r21 ).^2 ,3))) );
        ds12 =  nanmean( squeeze(sqrt(nansum( (r11 - r22).^2 ,3))) );
        ds21 =  nanmean( squeeze(sqrt(nansum( (r12 - r21).^2 ,3))) );
        ds22 =  nanmean( squeeze(sqrt(nansum( (r12 - r22).^2 ,3))) );
        %  we'll also want to swap the confidence scores
        conf11 = ref1_conf(1,tw,:); 
        conf12 = ref1_conf(2,tw,:); 
        % update ref1
        if ds12 + ds21 < ds11 + ds22 % ds12 < ds11 && ds21 < ds22
            ref1(1,tw,:) = r12;
            ref1(2,tw,:) = r11;
            ref1_conf(1,tw,:) = conf12; 
            ref1_conf(2,tw,:) = conf11; 
            if pars.verbose
                disp('swapped chn1 segments at ')
                disp([tw(1),tw(end)]);
            end
        end
    end
end

if ~isempty(posSwaps2)
    % if swapping sister identity at the trace
    posSwaps2 = [0; posSwaps2; tObs];
    for p=1:length(posSwaps2)-1
        tw = posSwaps2(p)+1:posSwaps2(p+1);
        % compute distance
        r11 = ref1(1,tw,:);
        r21 = ref2(1,tw,:);
        r12 = ref1(2,tw,:);
        r22 = ref2(2,tw,:);
        ds11 =  nanmean( squeeze(sqrt(nansum( (r11 - r21 ).^2 ,3))) );
        ds12 =  nanmean( squeeze(sqrt(nansum( (r11 - r22).^2 ,3))) );
        ds21 =  nanmean( squeeze(sqrt(nansum( (r12 - r21).^2 ,3))) );
        ds22 =  nanmean( squeeze(sqrt(nansum( (r12 - r22).^2 ,3))) );
        %  we'll also want to swap the confidence scores
        conf21 = ref2_conf(1,tw,:); 
        conf22 = ref2_conf(2,tw,:);
        if ds12 + ds21 < ds11 + ds22 %ds12 < ds11 && ds21 < ds22
            ref2(1,tw,:) = r22;
            ref2(2,tw,:) = r21;
            ref2_conf(1,tw,:) = conf22; 
            ref2_conf(2,tw,:) = conf21; 
            disp('swapped chn2 segments at ')
            disp([tw(1),tw(end)]);
        end
    end
end

%% assign identities based on length of data



n11 = nansum((ref1_conf(1,:,1)));
n12 = nansum((ref1_conf(2,:,1)));
n21 = nansum((ref2_conf(1,:,1)));
n22 = nansum((ref2_conf(2,:,1)));

% primary sisters are ref1(1,:,:) and ref2(1,:,:) chn1, chn2 respectively
% each channel has secondary sisters ref1(2,:,:) and ref2(2,:,:)
%
% Let's swap the primary/secondary labels so that the primary always have
% the longest traces. 

% ref1 and ref2 are now aligned sister pairs
if (n11 + n21) < (n12 + n22)
    ref1a = ref1;
    ref2a = ref2;
    ref1(1,:,:) = ref1a(2,:,:);
    ref1(2,:,:) = ref1a(1,:,:);
    ref2(1,:,:) = ref2a(2,:,:);
    ref2(2,:,:) = ref2a(1,:,:);
end


% %% show results
% cmapz = .5*hsv(zDepth);
% figure(4); clf;
% ts = 1:tObs;
% ax = [];
% alph = 0.2;
% ref_linewidth = 12; % line width
% dat_makersize = 2;
% dim = {['s=',num2str(s),'  x trace'],['s=',num2str(s),'  y trace']};
% for d=1:2
% ax(d) = subplot(2,1,d); 
%     % overlay the data
%     for z=1:zDepth
%         plot(ts(z:zDepth:end),squeeze(spotTrace1(s,z:zDepth:end,d,1)),'+','color',cmapz(z,:),'MarkerSize',dat_makersize); hold on;
%         plot(ts(z:zDepth:end),squeeze(spotTrace1(s,z:zDepth:end,d,2)),'x','color',cmapz(z,:),'MarkerSize',dat_makersize);
%         plot(ts(z:zDepth:end),squeeze(spotTrace2(s,z:zDepth:end,d,1)),'o','color',cmapz(z,:),'MarkerSize',dat_makersize);
%         plot(ts(z:zDepth:end),squeeze(spotTrace2(s,z:zDepth:end,d,2)),'s','color',cmapz(z,:),'MarkerSize',dat_makersize);
%     end
%     plot(ts,ref1(1,:,d),'-','LineWidth',ref_linewidth,'color',[0 0 1 alph]); hold on;
%     plot(ts,ref1(2,:,d),'-','LineWidth',ref_linewidth,'color',[0 1 1 alph]);
%     plot(ts,ref2(1,:,d),'-','LineWidth',ref_linewidth,'color',[1 0 0 alph]);
%     plot(ts,ref2(2,:,d),'-','LineWidth',ref_linewidth,'color',[1 1 0 alph]);
%     title(dim{d})
% end
% linkaxes(ax,'x');
% set(gcf,'color','w');


%% Use average traces to assign data
% a simple nearest ref should work at this point. 
%  we also add a max distance  


% minSister = 2.5;


ref1o = ref1; % preserve the original extrapolated traces
ref2o = ref2;

dis1 = sqrt(sum( (ref1(1,:,1:2) - ref1(2,:,1:2)).^2, 3)) ;
ref1(2,dis1<minSister,:) = inf;
dis2 = sqrt(sum( (ref2(1,:,1:2) - ref2(2,:,1:2)).^2, 3)) ;
ref2(2,dis2<minSister,:) = inf;

% figure(11); clf; plot(dis1,'.-'); hold on; plot(dis2,'.-'); legend('dis1','dis2')

% ref1 is chn1 which works for spotTrace1 (chn1)
dis_s1_rs1 = sqrt(sum((spotTrace1_s(:,1:2,1) - squeeze(ref1(1,:,1:2))  ).^2,2)); % sis 1
dis_s1_rs2 = sqrt(sum((spotTrace1_s(:,1:2,1) - squeeze(ref1(2,:,1:2))  ).^2,2)); % sis 2
dis_s2_rs1 = sqrt(sum((spotTrace1_s(:,1:2,2) - squeeze(ref1(1,:,1:2))  ).^2,2)); % sis 1
dis_s2_rs2 = sqrt(sum((spotTrace1_s(:,1:2,2) - squeeze(ref1(2,:,1:2))  ).^2,2)); % sis 2

spotTrace1_out_s = nan(size(spotTrace1_s));
is_sis1 = (dis_s1_rs1 <= dis_s1_rs2) & dis_s1_rs1 < maxDistToRef;
is_sis2 = (dis_s1_rs2 <= dis_s1_rs1) & dis_s1_rs2 < maxDistToRef; % need to avoid nans, can't use ~sis1
spotTrace1_out_s(is_sis1,:,1) = spotTrace1_s(is_sis1,:,1);  % anything in 1 that is marked as sis1
spotTrace1_out_s(is_sis2,:,2) = spotTrace1_s(is_sis2,:,1);  % anything in 1 that is marked as sis2


is2_sis1 = (dis_s2_rs1 <= dis_s2_rs2) &  dis_s2_rs1  < maxDistToRef; 
is2_sis2 = (dis_s2_rs2 <= dis_s2_rs1) &  dis_s2_rs2  < maxDistToRef; 
spotTrace1_out_s(is2_sis1,:,1) = spotTrace1_s(is2_sis1,:,2);  % anything in 2 that is marked as sis1
spotTrace1_out_s(is2_sis2,:,2) = spotTrace1_s(is2_sis2,:,2);  % anything in 2 that is marked as sis1


% ref2 is chn1 which works for spotTrace2 (chn2)
dis_s1_rs1 = sqrt(sum((spotTrace2_s(:,1:2,1) - squeeze(ref2(1,:,1:2)) ).^2,2)); % sis 1
dis_s1_rs2 = sqrt(sum((spotTrace2_s(:,1:2,1) - squeeze(ref2(2,:,1:2)) ).^2,2)); % sis 2
dis_s2_rs1 = sqrt(sum((spotTrace2_s(:,1:2,2) - squeeze(ref2(1,:,1:2)) ).^2,2)); % sis 1
dis_s2_rs2 = sqrt(sum((spotTrace2_s(:,1:2,2) - squeeze(ref2(2,:,1:2)) ).^2,2)); % sis 2

% figure(11); clf;
% plot(dis_s1_rs1,'r.'); hold on;
% plot(dis_s1_rs2,'b.'); legend('dis_s1_rs1','dis_s1_rs2')

spotTrace2_out_s = nan(size(spotTrace2_s));
is_sis1 = (dis_s1_rs1 <= dis_s1_rs2) & dis_s1_rs1 < maxDistToRef; 
is_sis2 = (dis_s1_rs2 <= dis_s1_rs1) & dis_s1_rs2 < maxDistToRef;   
spotTrace2_out_s(is_sis1,:,1) = spotTrace2_s(is_sis1,:,1);  % anything in 1 that is marked as sis1
spotTrace2_out_s(is_sis2,:,2) = spotTrace2_s(is_sis2,:,1);  % anything in 1 that is marked as sis2
is2_sis1 = (dis_s2_rs1 <= dis_s2_rs2) & dis_s2_rs1 < maxDistToRef;  
is2_sis2 = (dis_s2_rs2 <= dis_s2_rs1) & dis_s2_rs2 < maxDistToRef;  
spotTrace2_out_s(is2_sis1,:,1) = spotTrace2_s(is2_sis1,:,2);  % anything in 2 that is marked as sis1
spotTrace2_out_s(is2_sis2,:,2) = spotTrace2_s(is2_sis2,:,2);  % anything in 2 that is marked as sis1


% I think we should keep track of distance from reference as a quality
% score
%  In plenty of these traces (particularly those that borrow data from the
%  'sister' trace) we have a bunch of noise data far from reference that we
%  want to toss out.  


%% show results
if pars.figSisLabel > 0
    cmapz = .5*jet(zDepth);
    figSisLabel = figure(pars.figSisLabel); clf;
    figSisLabel.Position = [0 0 1200 1200];
    ts = 1:tObs;
    ax = zeros(1,4);
    alph = 0.2;
    ref_linewidth = 12; % line width
    dat_makersize = 4;
    dim = {['s=',num2str(s),'  x trace'],['s=',num2str(s),'  y trace']};
    k=0;
    for d=1:2
        k=k+1;
    ax(k) = subplot(4,1,k);
        plot(ts,ref1o(1,:,d),'-','LineWidth',ref_linewidth,'color',[0 0 1 alph]); hold on;
        plot(ts,ref1o(2,:,d),'-','LineWidth',ref_linewidth,'color',[0 1 1 alph]);
        plot(ts,ref2o(1,:,d),'-','LineWidth',ref_linewidth,'color',[1 0 0 alph]);
        plot(ts,ref2o(2,:,d),'-','LineWidth',ref_linewidth,'color',[1 1 0 alph]);
        % overlay the data
        for z=1:zDepth
            plot(ts(z:zDepth:end),squeeze(spotTrace1_s(z:zDepth:end,d,1)),'+','color',cmapz(z,:),'MarkerSize',dat_makersize); hold on;
            plot(ts(z:zDepth:end),squeeze(spotTrace1_s(z:zDepth:end,d,2)),'x','color',cmapz(z,:),'MarkerSize',dat_makersize);
            plot(ts(z:zDepth:end),squeeze(spotTrace2_s(z:zDepth:end,d,1)),'o','color',cmapz(z,:),'MarkerSize',dat_makersize);
            plot(ts(z:zDepth:end),squeeze(spotTrace2_s(z:zDepth:end,d,2)),'s','color',cmapz(z,:),'MarkerSize',dat_makersize);
        end
        plot(ts,squeeze(spotTrace1_s(:,d,1)),'b.','MarkerSize',dat_makersize); hold on;
        plot(ts,squeeze(spotTrace1_s(:,d,2)),'c.','MarkerSize',dat_makersize); hold on;
        plot(ts,squeeze(spotTrace2_s(:,d,1)),'r.','MarkerSize',dat_makersize); hold on;
        plot(ts,squeeze(spotTrace2_s(:,d,2)),'m.','MarkerSize',dat_makersize); hold on;
        colormap(cmapz); colorbar; clim([.5,zDepth+.5]); title(dim{d})
        title(['input ',dim{d}])
    end
    
    for d=1:2
        k=k+1;
    ax(k) = subplot(4,1,k);
    
        plot(ts,ref1o(1,:,d),'-','LineWidth',ref_linewidth,'color',[0 0 1 alph]); hold on;
        plot(ts,ref1o(2,:,d),'-','LineWidth',ref_linewidth,'color',[0 1 1 alph]);
        plot(ts,ref2o(1,:,d),'-','LineWidth',ref_linewidth,'color',[1 0 0 alph]);
        plot(ts,ref2o(2,:,d),'-','LineWidth',ref_linewidth,'color',[1 1 0 alph]);
    
        % overlay the data
        for z=1:zDepth
            plot(ts(z:zDepth:end),squeeze(spotTrace1_out_s(z:zDepth:end,d,1)),'+','color',cmapz(z,:),'MarkerSize',dat_makersize); hold on;
            plot(ts(z:zDepth:end),squeeze(spotTrace1_out_s(z:zDepth:end,d,2)),'x','color',cmapz(z,:),'MarkerSize',dat_makersize);
            plot(ts(z:zDepth:end),squeeze(spotTrace2_out_s(z:zDepth:end,d,1)),'o','color',cmapz(z,:),'MarkerSize',dat_makersize);
            plot(ts(z:zDepth:end),squeeze(spotTrace2_out_s(z:zDepth:end,d,2)),'s','color',cmapz(z,:),'MarkerSize',dat_makersize);
        end
        plot(ts,squeeze(spotTrace1_out_s(:,d,1)),'b.','MarkerSize',dat_makersize); hold on;
        plot(ts,squeeze(spotTrace1_out_s(:,d,2)),'c.','MarkerSize',dat_makersize); hold on;
        plot(ts,squeeze(spotTrace2_out_s(:,d,1)),'r.','MarkerSize',dat_makersize); hold on;
        plot(ts,squeeze(spotTrace2_out_s(:,d,2)),'m.','MarkerSize',dat_makersize); hold on;
        colormap(cmapz); colorbar; clim([.5,zDepth+.5]); title(dim{d})
        title(['final ',dim{d}])
    end 
    linkaxes(ax,'x');
    set(gcf,'color','w');
    if pars.saveFig
        SaveFigure(figSisLabel,'name',['s',num2str(s,'%03d'),'_SisLabel'],'formats',{'png'},'overwrite',true);
    end
end

%%
% figure(6); clf;
% for d=1:2
%     subplot(2,1,d);
%     alph = 0.2; ref_linewidth = 12;
%     plot(ts,ref1o(1,:,d),'-','LineWidth',ref_linewidth,'color',[0 0 1 alph]); hold on;
%     plot(ts,ref1o(2,:,d),'-','LineWidth',ref_linewidth,'color',[0 1 1 alph]);
%     plot(ts,ref2o(1,:,d),'-','LineWidth',ref_linewidth,'color',[1 0 0 alph]);
%     plot(ts,ref2o(2,:,d),'-','LineWidth',ref_linewidth,'color',[1 1 0 alph]);
% 
%     alph = 1; ref_linewidth = 1;
%     plot(ts,ref1(1,:,d),'.-','LineWidth',ref_linewidth,'color',[0 0 1 alph]); hold on;
%     plot(ts,ref1(2,:,d),'.-','LineWidth',ref_linewidth,'color',[0 1 1 alph]);
%     plot(ts,ref2(1,:,d),'.-','LineWidth',ref_linewidth,'color',[1 0 0 alph]);
%     plot(ts,ref2(2,:,d),'.-','LineWidth',ref_linewidth,'color',[1 1 0 alph]);
% 
% end