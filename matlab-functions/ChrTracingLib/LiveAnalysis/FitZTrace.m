function [zMat,stackHeight] = FitZTrace(pntArray,offTable,varargin)
% see also DaoFitZ
% Called by ProcessTraceGUI (dedicated?)
% 
% inputs 
%   pntArray matrix, N-spots x T-time x 6-datas: x,y,z,t,h,b 
% 
% outputs
%  fitsO - a structure with fields x,y,z,xm,ym,h,zw,za
%  many of the original x,y spots are now linked together in a z-stack 
%  and we combine them into a single x,y coordinates 
% 
% (x,y,z) - fitting centroid, x,y is from averaging all the x,y (2D) fits
% that went into the z
%  xm - x poistion from the linked 2D spot with the brightest fit
%  ym - y poistion from the linked 2D spot with the brightest fit
%  za - height of the z-gaussian above background, can be used to filter
%  weak spots
%  zw - sigma for the z-guassian - can be used to filter doublets etc.
% 


behavior = 'new';


if strcmp(behavior,'new')

defaults = cell(0,3);
defaults(end+1,:) = {'nmPerStep','positive',1000};  % 
defaults(end+1,:) = {'showPlots','boolean',true};  % 
defaults(end+1,:) = {'stackHeight','integer',0};  % 
defaults(end+1,:) = {'cycleStart','integer',0};  % 
defaults(end+1,:) = {'verbose','boolean',true};  %
defaults(end+1,:) = {'veryverbose','boolean',false};  %
pars = ParseVariableArguments(varargin,defaults,mfilename);


%% processs offset table to nm

% first, determine the z-scan cycle;
stageZ = offTable.offset;
if pars.stackHeight == 0 
    currZ = 0; oldZ = -inf; 
    [~,z0] =min(stageZ(1:100));
    z = z0;
    while currZ > oldZ
        oldZ = stageZ(z);
        z=z+1;
        currZ = stageZ(z);
    end
    stackHeight = z-z0;
    cycleStart = rem(z0,stackHeight);
else
    cycleStart = 1;
    stackHeight = pars.stackHeight;
end


figure(1); clf;
subplot(2,1,1); plot(stageZ(1:100),'.-');
title(['cycleStart = frame ',num2str(cycleStart),'.  stackHeight = ',num2str(stackHeight)'])
subplot(2,1,2); plot(stageZ,'.-');


% need to convert to nm.  1000 nm per step 
xx = stageZ(cycleStart:cycleStart+stackHeight-1); % just use the first 3 offset values, 
yy = 0:pars.nmPerStep:pars.nmPerStep*(stackHeight-1); %  pars.zscan_nm;% [-1000,0,1000]; % nm  (should become a parameter)
pp = polyfit(xx,yy,1);
yf = pp(1)*xx+pp(2);
stageZ_nm = pp(1)*stageZ+pp(2);  % offset-units to nm converted trajectory

% clean up stageZ
dZ = (diff(stageZ_nm)); % assumes a decreasing z value then a jump.
dropZ = dZ>2.5*pars.nmPerStep;
stageZ_nm([false; dropZ]) = stageZ_nm([dropZ; false])+pars.nmPerStep;

figure(1); clf;
subplot(2,1,1); plot(stageZ_nm(1:100),'.-');
title(['cycleStart = frame ',num2str(cycleStart),'.  stackHeight = ',num2str(stackHeight)'])
subplot(2,1,2); plot(stageZ_nm,'.-');



%%
    % Gaussian fit the brightness and read that in as z. 
    ftype = fittype('exp(-((x-mu_x)/(2*500)).^2)*a',...
                     'coeff', {'a','mu_x'},'ind',{'x'});
    [nPts,tTotal,~] = size(pntArray);
    tSteps = tTotal/stackHeight;
    zMat = nan(nPts,tTotal);
    k=0; % counter for the total number of fits computed
    for n=1:nPts  % n=2;  zMat = nan(nPts,tTotal);
         if pars.showPlots
            figure(2); clf;
            ax1 = subplot(2,1,1); plot(stageZ_nm,'.-');  title(n);
            ylabel('stage-Z (nm)');
            ax2 = subplot(2,1,2); 
            plot(squeeze(pntArray(n,:,5)),'.-'); hold on;
            plot(squeeze(pntArray(n,:,6)),'.-'); hold on;
            ylabel('brightness'); xlabel('time');
            linkaxes([ax1,ax2],'x');
         end

        if pars.verbose
            if rem(n,2)==0
                disp(['fitting ',num2str(n/nPts*100,3),'% complete'])
            end
        end
        bkd = nanmedian(pntArray(n,:,6));
        for t=1:tSteps
            t1 = stackHeight*(t-1)+1;
            t2 = stackHeight*t;
            b = pntArray(n,t1:t2,5)'  ; % + pntArray(n,t1:t2,6)'; % height plus background  
            nObs = stackHeight - sum(isnan(b));
            if nObs == 0 % skip if there is no data
                continue
            elseif nObs == 1
                % zMat(n,t1:t2) = stageZ_nm(t1-1+~isnan(b));
                continue
            end
            b(isnan(b)) = 0;
            z=stageZ_nm(t1:t2);  % we have some data to fit
            % squeeze(pntArray(n,t1:t2,1:6))

             try
               %    Gaussian fit
               [~,i] = max(b);
               sp = [max(b),    z(i)];
               lb = [   0,         0]; 
               ub = [4*max(b),  max(z)];
        
                fitZ = fit(z,b,ftype,...
                     'StartPoint',sp,...
                     'Lower',lb,...
                     'Upper',ub,...
                     'TolFun',1e-6,'TolX',1e-6,'MaxIter',4000,'MaxFunEvals',600);
               %  figure(10); clf; plot(fitZ,z,b); pause(.1); xlabel('z'); ylabel('brightness'); 
                zMat(n,t1:t2) = fitZ.mu_x;
                k=k+1;

                ci = confint(fitZ,.75);
                if (ci(2,2)-ci(1,2))/(max(z)-min(z)) > 1
                    poorSpot = ['detected poor fit on spot = ',num2str(n),' t=',num2str(t)];
                    if pars.veryverbose
                        disp(poorSpot);
                    end
                    zMat(n,t1:t2) = nan;
                    % figure(10); clf; plot(fitZ,z,b); pause(.1); xlabel('z'); ylabel('brightness'); title(poorSpot);
                    % disp('place debug here')
                end
                if  rem(k,1000) == 0 && pars.showPlots
                  figure(3); clf; plot(fitZ,z,b); 
                  pause(.1); xlabel('z (nm)'); ylabel('brightness'); 
                  title(['example fit on spot = ',num2str(n),' t=',num2str(t)]);
                  pause(.01); 
                end

            catch er
                disp(er.getReport);
                disp('place debug here')
             end   
        end
        if pars.showPlots
            figure(1); clf; plot(zMat(n,:),'.-');   title(n); 
            xlabel('time'); ylabel('z (nm)'); pause(.01);  
        end
    end






























else
%% Old version 
% we use z-position measured with the IR laser for the frame and the
% z-stack data.  Though it's impressive there are so many paired non-zero
% detection events with a 1 um step size. 

% z-fitting
defaults = cell(0,3);
defaults(end+1,:) = {'zFrames','integer',3}; % 
defaults(end+1,:) = {'zscan_nm','array',[-1000,0,1000]}; % 
defaults(end+1,:) = {'smoothZ','boolean',false}; % use smoothing?

pars = ParseVariableArguments(varargin,defaults,mfilename);


% offFile =  regexprep(daxName1,'_C1.dax','.off');
% offTable = ReadTableFile(offFile,'delimiter',' ');
    
    [nTraces,nFrames,~] = size(pntArray1a);
    
    
    
    o = [offTable.offset(2:end); nan]; % shifted by 1 frame for some reason
    o(o==0) = nan;
    o(o>quantile(o,.9999)) = nan;
    % figure(3); clf; plot(o(1:10),'.-');  % quick look at the offset oscillations
    
    % need to convert to nm.  1000 nm per step 
    % xx = [quantile(o,.15), quantile(o,.5), quantile(o,.85)]; % offset units
    xx = o(1:pars.zFrames); % just use the first 3 offset values, 
    % xx = [median(o(1:pars.zFrames:end)),median(o(2:pars.zFrames:end)),median(o(3:pars.zFrames:end))]  
    yy =  pars.zscan_nm;% [-1000,0,1000]; % nm  (should become a parameter)
    pp = polyfit(xx,yy,1);
    yf = pp(1)*xx+pp(2);
    onm = pp(1)*o+pp(2);  % offset-units to nm converted trajectory
    
    % figure(1); clf; plot(xx,yy,'.');
    % hold on; plot(xx,yf,'-');

zMat = nan(nTraces,nFrames);
for s = 1:nTraces 
    try
    x = squeeze(pntArray1a(s,:,1));
    y = squeeze(pntArray1a(s,:,2));
    h = squeeze(pntArray1a(s,:,5));
    b = squeeze(pntArray1a(s,:,6));
    
    % clean up background
    bkd = b;  
    % leading and trailing nans screw up the linear interpolation. we need to elinate these first.
    firstNonNan = find(~isnan(bkd),100,'first');
    lastNonNan = find(~isnan(bkd),100,'last');
    bkd(1) = mean(bkd(firstNonNan));
    bkd(end) = mean(bkd(lastNonNan));%
    bkd = fillmissing(bkd,'linear');
    bkd = smooth(bkd,.2); % photobleaching, gradual decay
    % figure(1); clf; plot(bkd); ylabel('bkd'); title('background')
    
    zDat = h;
    zDat(isnan(h)) = bkd(isnan(h));
    
    % figure(1); clf; 
    % ax1 = subplot(2,1,1); plot(zDat(fs),'.-'); hold on; plot(bkd(fs),'--');
    % ax2 = subplot(2,1,2); plot(o(fs),'.-');
    % linkaxes([ax1,ax2],'x')
    
    
    bkdMin = quantile(bkd,.001);
    bkdMax = quantile(bkd,.999);
    sigMin = quantile(h,.001) - bkdMin;
    sigMax = quantile(h,.999) + bkdMax;
    zf = nan(1,nFrames);
     for t=2:pars.zFrames:nFrames-2
         xx = onm(t:t+2);
        zDat_t = zDat(t:t+2)';
        if sum(zDat_t == bkd(t:t+2))<=1  % no missed data or 1 missed data (3 obs)
            [v,i] = min(zDat(t:t+2)); % just use 2 brightest
            hb= zDat_t./bkd(t:t+2) -1;     
            if i==1 % ignore first position
                zf(t) = hb(2)/sum(hb(2:3))*xx(2) + hb(3)/sum(hb(2:3))*xx(3);
            elseif i==3 % ignore last position
                zf(t) = hb(1)/sum(hb(1:2))*xx(1) +  hb(2)/sum(hb(1:2))*xx(2)  ;
            else % is middle is missed and outsides are detected, something is wrong,  
                zf(t) =  nan;
            end
        elseif sum(zDat_t == bkd(t:t+2))<1  % only 1 is detected, best guess is that
            zf(t) = xx(t); 
        end
        % if rem((t),100)==0
        %     disp([num2str(t/nFrames*100),'% complete']); pause(.01);
        % end
     end
    
    % figure(4); clf; 
    % subplot(3,1,1); plot(zDat(fs),'.-'); hold on; plot(bkd(fs)); ylabel('z bright')
    % subplot(3,1,2); plot(o(3:4:end),'.-'); ylabel('offset')
    % subplot(3,1,3); plot(zf(fs),'.-'); ylabel('z-fit');%  ylim([-3 3])
    
    zfMax = quantile(zf,.95);
    zfMin = quantile(zf,.05);
    zff = zf; 
    zff(zf==0) = nan;
    zff(zf>zfMax | zf<zfMin) = nan;
    if pars.smoothZ
        try
            zs = smooth(zff,5,'rlowess');  % intervals 
        catch er
            warning(er.message);
            warning(er.getReport);
            zs = zff;
        end
    else
        zs = zff;
    end
    zs(isnan(x)) = nan;

    zMat(s,:) = zs;
    disp([num2str(s/nTraces*100),'% complete']); pause(.01);
    catch er
        warning(er.getReport);
        
    end

end


% figure(2); clf; 
% a1 = subplot(3,1,1); plot((x(fs)-nanmean(x(fs)))*npp,'.-'); ylabel('x (nm)')
% a2 = subplot(3,1,2); plot((y(fs)-nanmean(y(fs)))*npp,'.-'); ylabel('y (nm)')
% a3 = subplot(3,1,3); 
% plot(zs,'.-'); yl = get(gca,'ylim');
% hold on; plot(zff(fs),'.'); ylabel('z-fit (nm)');%  ylim([-3 3])
% ylim(yl); title('z-fit and z-interp.')
% linkaxes([a1,a2,a3],'x');

end

