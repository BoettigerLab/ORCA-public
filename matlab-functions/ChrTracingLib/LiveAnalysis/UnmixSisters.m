function [spotTrace1_out,spotTrace2_out] = UnmixSisters(spotTrace1,spotTrace2,varargin)
%% Approach
% 1. look for jumps in the relative distance between each sister on trace 1
% compared to a un-gapped reference version of trace 2. Then swap.
% The *relative* measure is important because the cells  /chromosomes move a
% lot more than the relative distances between spots on the same chromosome.
% Any rough approximation of relative distance will be an improvment (and 
% we apply the identical comparison to both sisters, so there's no
% distortions).  Because the other traces may not be seen in the same
% z-planes, and may contain small gaps, we use an interpolated version just
% as handle into the chromosome reference frame.  This interpolation is not
% propogated and doesn't need to be perfect. 
% 2. before accepting the swap, we make sure the sister trace also looks
% like real data and not just noise.  (because we literaly just save the
% first point in the FOV not in the trace, 
%
% Things still to address in subsequent functions
%    - no guarentee that different z all select the same sister
%    
%% Notes
%  3/29/24 - still in development 

%%
defaults = cell(0,3);
defaults(end+1,:) = {'figSis','integer',4}; % 0 for off
defaults(end+1,:) = {'ref_window','positive',40}; % should be at least size z-depth. Better to be 4 or 5x z-depth
defaults(end+1,:) = {'maxJmp','positive',10}; %  maximum spot-to-spot detection jump, relative to spot-to-spot step, to be considered an outlier and a potential sister swap.  this effects how often we try to swap sisters
defaults(end+1,:) = {'w','positive',1000}; % window (in total frames) to compare pre and post distances over.  
defaults(end+1,:) = {'verbose','boolean',false};
pars = ParseVariableArguments(varargin,defaults,mfilename);

pars.figSis = 4;
ref_window = 40; % should be at least size z-depth. Better to be 4 or 5x z-depth
% This effects how much 
% % interpolation we are willing to do to estimate the position
% of trace 1 (sister i) relative to trace 2. 
maxJmp = 10; % this effects how often we try to swap sisters
w = 1000; % how far foward and back to look for trajectory

minSisLength = 0.33; % at least x fraction the length of the sister 
maxSisStepVar = 1.5;
%%
spotTrace1_out = spotTrace1;
spotTrace2_out = spotTrace2; 

[nSpots,tObs,dDims] = size(spotTrace1_out);
zDepth = pars.zDepth;
ts = 1:tObs;
idx = 1:(tObs/zDepth); %#ok<NASGU>  (used for plotting)


% could split these into cell arrays first and par for this if needed
for s=1:nSpots % s=132    
         % in case the 2nd color is in a different z-plane, we interpolate first
         %    approach 1, interpolating all z, creates problems when the
         %    reference has 2 sisters that are in different z as well as
         %    have different x,y. 

        % tr2_best = nan(zDepth,tObs/zDepth);
        % for z=1:zDepth
        % tr2_best(z,:) = fillmissing(squeeze(spotTrace2_out(s,z:zDepth:end,1,1)),'movmean',ref_window,'maxGap',ref_window);
        % end
        % [nTr2,z2] = max(sum(~isnan(tr2_best),2));
        % x2 = fillmissing(squeeze(spotTrace2_out(s,z2:zDepth:end,1,1)),'movmean',ref_window,'maxGap',ref_window);
        % y2 = fillmissing(squeeze(spotTrace2_out(s,z2:zDepth:end,2,1)),'movmean',ref_window,'maxGap',ref_window);
        % xy2 = [x2;y2]';


    currStack2 = nan(nFrames/zDepth,8,zDepth);
    for z =1:zDepth
        currStack2(:,:,z) = squeeze(spotTrace2(s,z:zDepth:end,:,1));
    end
    ObsPerZ_2 = squeeze(sum(~isnan(currStack2(:,d,:)),1));
    ref_trace = max( currStack2(:,1:2,:),[],3); % max ensures that we only take 1 from each z-trace, and that we tend to take a consitent one   
    ref_trace(:,1) = fillmissing(ref_trace(:,1) ,'movmean',ref_window,'maxGap',ref_window);
    ref_trace(:,2) = fillmissing(ref_trace(:,2) ,'movmean',ref_window,'maxGap',ref_window);
    xy2 = ref_trace;


    for z = 1:zDepth % loop over z  z =3
        showSwap = false;
        % First we compute the jumps in relative distance of each sister in trace 1
        % from its Trace2 spot. To ensure alignment, we use an interpolated version
        % of trace2.  
        tz = ts(z:zDepth:end); % keep track of actual frame
        tt = tz;
        x1 = squeeze(spotTrace1_out(s,z:zDepth:end,1,1));
        y1 = squeeze(spotTrace1_out(s,z:zDepth:end,2,1));
        
        
        x1r = x1 - x2; % relative x
        y1r = y1 - y2; % rleative y

        gaps = isnan(x1r);
        x1r(gaps) = [];
        y1r(gaps) = [];
        tt(gaps) = [];
        
        jmp1 = sqrt( (x1r(1:end-1) - x1r(2:end)).^2 +  (y1r(1:end-1) - y1r(2:end)).^2 );
        jmp5 = sqrt( (x1r(1:end-5) - x1r(6:end)).^2 +  (y1r(1:end-5) - y1r(6:end)).^2 );
        
        % % just for troubleshooting
        % figure(4); clf; 
        % plot(tt(1:end-1),jmp1,'o'); hold on;
        % hold on; plot(tt(1:end-5),jmp5,'o');
        % 
        % figure(4); clf;
        % yyaxis left;
        % plot(ts(z:zDepth:end),spotTrace1_out(s,z:zDepth:end,1,1),'b.'); hold on;
        % plot(ts(z:zDepth:end),spotTrace1_out(s,z:zDepth:end,1,2),'r.'); hold on;
        % yyaxis right;
        % plot(tt(1:end-1),jmp1,'o'); hold on;
        % plot(tt(2:end),jmp1,'o'); hold on;
        
        tt_up = tt(1:end-1);
        tt_down = tt(2:end);
        t_jmp = [tt_up(jmp1 > maxJmp*median(jmp1_up)); tt_down(jmp1 > maxJmp*median(jmp1_up))];
        t_jump1z = 1+floor(t_jmp/zDepth);
        nJumps = size(t_jump1z,2);
        
        % if pars.figSis > 0
        %     figure(pars.figSis); clf;
        %     yyaxis left;
        %     plot(tz,spotTrace2(s,z2:zDepth:end,1,1),'.','color',.7*ones(1,3)); hold on;
        %     plot(tz,spotTrace1_out(s,z:zDepth:end,1,1),'b.'); hold on;
        %     plot(tz,spotTrace1_out(s,z:zDepth:end,1,2),'r.'); hold on;
        %     yyaxis right;
        %     plot(tz(t_jump1z),ones(length(t_jump1z),1),'ko'); hold on;
        %     legend('trace2',['trace1 sis1, z=',num2str(z)],['trace1 sis2, z=',num2str(z)],'detected jumps','location','Best')
        % end
    
        xy_sis1 = squeeze(spotTrace1_out(s,z:zDepth:end,1:2,1)) - xy2;
        xy_sis2 = squeeze(spotTrace1_out(s,z:zDepth:end,1:2,2)) - xy2;
    
        % before we get into the computationally slow stuff, let's check
        % that we have enough sister data in the trace to be worth
        % processing
        len1 = length(x1r);
        len2 = sum(~isnan(xy_sis2(:,1)));
        if len2 > minSisLength*len1
            % clean-up the sister trace (This is slow!)
            xy_sis2_filt = RemoveJumps(xy_sis2,'localRegion',10,'maxDiff',6,'maxAbsStep',3,'removeLoners',false);
            % % for troubleshooting, view the sister trace jumps
            % figure(4); clf; subplot(2,1,1);
            %  plot(xy_sis2(:,2),'o'); hold on;
            %  plot(xy_sis2_filt(:,2),'.-');
            %  subplot(2,1,2); plot(xy2(:,1),'.-');
            % after plotting testing, we accept the cleaned up version
            xy_sis2 = xy_sis2_filt; 
        
            
            % test if sis2 looks real
            gaps2 = isnan(xy_sis2(:,1));
            tt_sis = tz;
            x_sis2_test = xy_sis2(:,1);
            x_sis2_test(gaps2,:) = [];
            tt_sis(gaps2) = [];
            stp1 = nanmedian(abs(diff(x1r)));
            stp2 = nanmedian(abs(diff(x_sis2_test)));
            len1 = length(x1r);
            len2 = length(x_sis2_test);

            % figure(3); clf;
            % plot(x1r,'.'); hold on; plot(x_sis2_test,'.');
            
            valid_sister = stp2 < maxSisStepVar*stp1  & len2 > minSisLength*len1 ;
            if valid_sister
                for t=1:nJumps 
                    tu0 = max(t_jump1z(1,t)-w,1); % upstream start
                    tu1 = t_jump1z(1,t); % upstream end
                    td0 = t_jump1z(2,t); % downstream start
                    td1 = min(t_jump1z(2,t)+w,tObs/zDepth); % downstream end
                    sis1_pre = xy_sis1(tu0:tu1,1:2);
                    sis1_post = xy_sis1(td0:td1,1:2);
    
                    % % just for troubleshooting
                    % %  figure(4); clf;
                    % sis2_pre = xy_sis2(tu0:tu1,1:2);
                    % sis2_post = xy_sis2(td0:td1,1:2);
                    % plot(idx(t_jump1z(1,t)-w:t_jump1z(1,t)), sis1_pre(:,1),'.'); hold on;
                    % plot(idx(t_jump1z(2,t):t_jump1z(2,t)+w), sis1_post(:,1),'.'); hold on;
                    % plot(idx(t_jump1z(1,t)-w:t_jump1z(1,t)), sis2_pre(:,1),'.'); hold on;
                    % plot(idx(t_jump1z(2,t):t_jump1z(2,t)+w), sis2_post(:,1),'.'); hold on;
                    % legend('sis1-pre','sis1-post','sis2-pre','sis2-post')
                
                    preJump_spotDist =  sqrt(sum(nanmedian(abs(sis1_pre )).^2));
                    postJump_spotDist = sqrt(sum(nanmedian(abs(sis1_post)).^2));
                
                    sis1_swap = xy_sis1;
                    sis1_swap(1:t_jump1z(1,t),:) = xy_sis1(1:t_jump1z(1,t),:);
                    sis1_swap(t_jump1z(2,t):end,:) = xy_sis2(t_jump1z(2,t):end,:);
                    sis2_swap = xy_sis2;
                    sis2_swap(1:t_jump1z(1,t),:) = xy_sis2(1:t_jump1z(1,t),:);
                    sis2_swap(t_jump1z(2,t):end,:) = xy_sis1(t_jump1z(2,t):end,:);
                    
                    swap1_pre = sis1_swap(tu0:tu1,1:2);
                    swap1_post = sis1_swap(td0:td1,1:2);
                    swap_preJump_spotDist =  sqrt(sum(nanmedian(abs(swap1_pre )).^2));
                    swap_postJump_spotDist = sqrt(sum(nanmedian(abs(swap1_post)).^2));
                
                    % If the distance changes less abruptly after switchnig the sister
                    % assignments at t_jump, we accept the switch. That is, if pre-jump
                    % they were close and swapping makes the post-jump Tet/CuO pair close,
                    % we accept it, if they were far and swapping keeps them far, we keep
                    % that - we assume no rapid change in 3D separation between pairs if
                    % changing sister identities would remove the change. 
                    if abs(postJump_spotDist-preJump_spotDist) > abs(swap_postJump_spotDist-swap_preJump_spotDist)
                        disp(['swapping sisters at t= ',num2str(t_jmp(:,t)' )] )
                        xy_sis1 = sis1_swap;
                        xy_sis2 = sis2_swap;
                        showSwap = true;
                    end
                end
                % insert the new data back into 'spotTrace'
                %    for simplicity, keep the sister that is physically
                %    closer to the reference as sis1 and other as sis2.
                %    Later steps will explore swapping
                if nanmean(abs(xy_sis1(:,1))) <= nanmean(abs(xy_sis2(:,1)))
                    spotTrace1_out(s,z:zDepth:end,1:2,1) = xy_sis1 + xy2;
                    spotTrace1_out(s,z:zDepth:end,1:2,2) = xy_sis2 + xy2;
                else
                    spotTrace1_out(s,z:zDepth:end,1:2,1) = xy_sis2 + xy2;
                    spotTrace1_out(s,z:zDepth:end,1:2,2) = xy_sis1 + xy2;
                end
            end
        
            if pars.figSis > 0  && showSwap
                figure(pars.figSis); clf; 
                subplot(2,2,1); d = 1;
                for zz2=1:zDepth; plot(tz, squeeze( spotTrace1(s,zz2:zDepth:end,d,1) ),'.-','color',[.8,.8,1]); hold on; end
                for zz2=1:zDepth; plot(tz, squeeze( spotTrace1(s,zz2:zDepth:end,d,2) ),'.-','color',[.8,1,1]); hold on; end
                for zz2=1:zDepth; plot(tz, squeeze( spotTrace2(s,zz2:zDepth:end,d,1) ),'.-','color',[1,.8,.8]); hold on; end
                for zz2=1:zDepth; plot(tz, squeeze( spotTrace2(s,zz2:zDepth:end,d,2) ),'.-','color',[.8,1,1]); hold on; end
                % p1 = plot(tz, squeeze( spotTrace2(s,z2:zDepth:end,d,1) ),'r.-'); hold on;
                p1 = plot(tz, squeeze( xy2(:,d) ),'r.-'); hold on;
                p2 = plot(tz, squeeze( spotTrace1(s,z:zDepth:end,d,1) ),'b.-'); hold on;
                p3 = plot(tz, squeeze( spotTrace1(s,z:zDepth:end,d,2) ),'c.-'); hold on;
                title(['s=',num2str(s), '  z=',num2str(z),'  orig x']);
                legend([p1,p2,p3],'trace2',['trace1 sis1, z=',num2str(z)],['trace1 sis2, z=',num2str(z)],'location','Best') 
                subplot(2,2,3); d = 2;% 'y
                for zz2=1:zDepth; plot(tz, squeeze( spotTrace1(s,zz2:zDepth:end,d,1) ),'.-','color',[.8,.8,1]); hold on; end
                for zz2=1:zDepth; plot(tz, squeeze( spotTrace2(s,zz2:zDepth:end,d,1) ),'.-','color',[1,.8,.8]); hold on; end
                % p1 = plot(tz, squeeze( spotTrace2(s,z2:zDepth:end,d,1) ),'r.-'); hold on;
                p1 = plot(tz, squeeze( xy2(:,d) ),'r.-'); hold on;
                p2 = plot(tz, squeeze( spotTrace1(s,z:zDepth:end,d,1) ),'b.-'); hold on;
                p3 = plot(tz, squeeze( spotTrace1(s,z:zDepth:end,d,2) ),'c.-'); hold on;
                title('orig y')
                subplot(2,2,2); d = 1;
                for zz2=1:zDepth; plot(tz, squeeze( spotTrace1_out(s,zz2:zDepth:end,d,1) ),'.-','color',[.8,.8,1]); hold on; end
                for zz2=1:zDepth; plot(tz, squeeze( spotTrace2(s,zz2:zDepth:end,d,1) ),'.-','color',[1,.8,.8]); hold on; end
                % p1 = plot(tz, squeeze( spotTrace2(s,z2:zDepth:end,d,1) ),'r.-'); hold on;
                p1 = plot(tz, squeeze( xy2(:,d) ),'r.-'); hold on;
                p2 = plot(tz, squeeze( spotTrace1_out(s,z:zDepth:end,d,1) ),'b.-'); hold on;
                p3 = plot(tz, squeeze( spotTrace1_out(s,z:zDepth:end,d,2) ),'c.-'); hold on;
                title(['s=',num2str(s), '  z=',num2str(z),'  new x']);
                legend([p1,p2,p3],'trace2',['trace1 sis1, z=',num2str(z)],['trace1 sis2, z=',num2str(z)],'location','Best') 
                subplot(2,2,4); d = 2;
                for zz2=1:zDepth; plot(tz, squeeze( spotTrace1_out(s,zz2:zDepth:end,d,1) ),'.-','color',[.8,.8,1]); hold on; end
                for zz2=1:zDepth; plot(tz, squeeze( spotTrace2(s,zz2:zDepth:end,d,1) ),'.-','color',[1,.8,.8]); hold on; end
                % p1 = plot(tz, squeeze( spotTrace2(s,z2:zDepth:end,d,1) ),'r.-'); hold on;
                p1 = plot(tz, squeeze( xy2(:,d) ),'r.-'); hold on;
                p2 = plot(tz, squeeze( spotTrace1_out(s,z:zDepth:end,d,1) ),'b.-'); hold on;
                p3 = plot(tz, squeeze( spotTrace1_out(s,z:zDepth:end,d,2) ),'c.-'); hold on;
                title('new y')
                set(gcf,'color','w');
                pause();
            end
        end
    end
end
%%
for s=1:nSpots
   figure(3); clf; 
        ax = [];
        zsym = {'<','^','>','v','s'};
        w=30;
        for d=1:2
            % ref1(:,d) = fillmissing(ref1(:,d),'linear','maxGap',30);
            % ref2(:,d) = fillmissing(ref2(:,d),'linear','maxGap',30);
            
            for z=1:zDepth
                ref1 = squeeze(spotTrace1_out(s,z:zDepth:end,1:3,1));
                ref2 = squeeze(spotTrace2_out(s,z:zDepth:end,1:3,1));
                ref1(:,d) = fillmissing(ref1(:,d),'movmean',w,'maxGap',w);
                ref2(:,d) = fillmissing(ref2(:,d),'movmean',w,'maxGap',w);
                ref1(:,d) = smooth(ref1(:,d),'moving',w);
                ref2(:,d) = smooth(ref2(:,d),'moving',w);
                figure(3); 
                ax(d) = subplot(2,1,d); plot(ts(z:zDepth:end),ref1(:,d),'r-');
                hold on; plot(ts(z:zDepth:end),ref2(:,d),'b-');
               

                plot(ts(z:zDepth:end),squeeze(spotTrace1_out(s,z:zDepth:end,d,1)),['r',zsym{z}]);
                plot(ts(z:zDepth:end),squeeze(spotTrace1_out(s,z:zDepth:end,d,2)),['m',zsym{z}]);
                plot(ts(z:zDepth:end),squeeze(spotTrace2_out(s,z:zDepth:end,d,1)),['b',zsym{z}]);
                plot(ts(z:zDepth:end),squeeze(spotTrace2_out(s,z:zDepth:end,d,2)),['c',zsym{z}]);
            end
        end
        linkaxes(ax,'x'); title(['s=',num2str(s)]);
        pause
end