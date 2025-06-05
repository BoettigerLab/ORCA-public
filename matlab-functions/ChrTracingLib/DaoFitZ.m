function fitsO = DaoFitZ(fits1,varargin)
% inputs 
%   fits1 - a structure from LoadHD5Fits.
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

defaults = cell(0,3);
defaults(end+1,:) = {'sigma','positive',5};  % width
defaults(end+1,:) = {'sigmaRange','positive',30}; % threshold in brightness units  
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'selFrames','float',0}; % select frames.  0 for all
defaults(end+1,:) = {'seedThresh','float',0.0}; % autoselect threshold for grouping localizations in z. [OBSOLETE]

defaults(end+1,:) = {'minSpotsPerTrace','positive',4}; % autoselect threshold for grouping localizations in z. 
defaults(end+1,:) = {'seedBinResolution','positive',2.5}; % autoselect threshold downsampling for grouping localizations in z. 
defaults(end+1,:) = {'showSeed','boolean',false}; % autoselect threshold downsampling for grouping localizations in z. 
defaults(end+1,:) = {'seedLoci','float',[]};  
defaults(end+1,:) = {'maxStep','positive',3}; % 
defaults(end+1,:) = {'maxDistToSeed','positive',6}; % 
defaults(end+1,:) = {'symbol','string','.'}; 
defaults(end+1,:) = {'z_per_frame','array',[]}; % 
defaults(end+1,:) = {'figHandle','freeType',[]}; % 
defaults(end+1,:) = {'multiFit','boolean',false}; % 
defaults(end+1,:) = {'multiFit_res','positive',2}; % resolution (in pixels) at which to attempt distinguishing
defaults(end+1,:) = {'multiFit_thetaAbs','fraction',200}; %  threshold in absolute brightness units  
defaults(end+1,:) = {'multiFit_thetaMax','positive',.75}; % threshold relative to max brightness along the z-pixel stack
defaults(end+1,:) = {'multiFit_zWindow','positive',1.4}; % threshold in brightness units  
defaults(end+1,:) = {'output',{'struct','table'},'struct'}; % update to table in future -- currently just struct for some backwards compatibility 
defaults(end+1,:) = {'showPoorFits','boolean',false}; %good for troubleshooting
pars = ParseVariableArguments(varargin,defaults,mfilename);

showPlot = false; % just for troubleshooting

% folderCal2 = 'M:\2023-08-26_130kb_rad21_dTag\';
% daxName = [folderCal2,'beads_0001_C1.dax'];
% [daxIm,daxInfo] = LoadDax(daxName);
% [h,w,zs] = size(daxIm);
% binName = regexprep(daxName,'.dax','.hdf5');
% fits1 = LoadHD5Fits(binName);

if istable(fits1)
     fitsO = DaoFitZ_FromTable(fits1,'parameters',pars)
else
    
    xs = {fits1.x};
    ys = {fits1.y};
    if isempty(pars.seedLoci)
    %     xy0 = TraceCenters(xs,ys,'autoSelectThreshold',pars.seedThresh,'binResolution',pars.seedBinResolution,'showPlots',pars.showSeed);
        xy0 = TraceCenters(xs,ys,'minSpotsPerBin',pars.minSpotsPerTrace,'minSpotsPerTrace',pars.minSpotsPerTrace, ...
            'binResolution',pars.seedBinResolution,'showPlots',pars.showSeed);
    else
        xy0 = pars.seedLoci;
    end

    
    %%
    if pars.selFrames == 0
        try
        selFrames = [fits1(:).frame]-[fits1(1).frame]+1; % must start from 1 as these are used to index out of the cell xs{}  
        catch
            maxFrame = max([fits1.frame]);
            minFrame = min([fits1.frame]);
            selFrames = 1:(maxFrame-minFrame+1);
        end
    else
        selFrames = pars.selFrames;
    end
    nFrames = length(selFrames);
    cmap = hsv(nFrames);
    
    tic
    nPts = size(xy0,1); % length(xs{1});
    idx = cell(nFrames,1); 
    pntChain = nan(nPts,2,nFrames+1);
    noChains = true;
    for f=1:nFrames % -1
        % disp(f)
        try
            fn = selFrames(f);
            if pars.verbose
                if rem(f,1000)==0
                    disp(['linking ',num2str(f/nFrames*100,3),'% complete'])
                end
            end
            if noChains    
                xy1 = xy0;   % seed point 
                xy2 = [xs{fn},ys{fn}];   % first data frame 
    
                if ~isempty(xy2)
                    [matched,cost] = MatchPoints(xy2,xy1,'maxDist',pars.maxDistToSeed); % s 
                    m =matched(cost<pars.maxDistToSeed,:);  
                else
                    m = [];
                end
                if ~isempty(m)
                    nMatched = size(m,1);
                    pntChain(m(:,2),1:2,f) = xy2(m(:,1),1:2); % indexed in order of starting point
                    noChains = false;
                else
                    noChains = true;
                end
            else
                xy1 = pntChain(:,1:2,f-1); % this is 1 frame ahead of the data, because we kept the seed
                xy2 = [xs{fn},ys{fn}];   % first data frame 
                
                % if missed in last frame, use previous frame
                missed = isnan(xy1(:,1));
                ff = f-2;
                while any(missed) && (ff >=1)
                    xy1(missed,:) = pntChain(missed,1:2,ff);
                    missed = isnan(xy1(:,1));
                    ff=ff-1;
                end
                xy1n = xy1;  % this is the new xy1, which hassing missing values back filled from the last-spot seen in the trace.  
                %  xy1n(isnan(xy1(:,1)),:) = 0; % if any are still missing, replace with 0
                xy1n(isnan(xy1(:,1)),:) = xy0(isnan(xy1(:,1)),:); % if still missing, replace with seed.  
                [matched,cost] = MatchPoints(xy2,xy1n,'maxDist',pars.maxStep);
                m =matched(cost<pars.maxStep,:);      
                if ~isempty(m) 
                    pntChain(m(:,2),1:2,f-1) = xy1(m(:,2),1:2);
                    pntChain(m(:,2),1:2,f) = xy2(m(:,1),1:2); % indexed in order of starting point
                end
            end
            idx{f} = m;
                % data(m(:,2),all_data_types,frameNum) = dataInFrame_f(m(:,1),all_data_types) 
            if showPlot && f>2
                figure(pars.figHandle);
                x1 = squeeze(pntChain(:,1,f-1))';
                x2 = squeeze(pntChain(:,1,f))';
                y1 = squeeze(pntChain(:,2,f-1))';
                y2 = squeeze(pntChain(:,2,f))';
                plot([x1;x2],[y1;y2],'color',cmap(f,:)); hold on;
                plot([x1;x2],[y1;y2],pars.symbol,'color',cmap(f,:));
                plot(xy2(:,1),xy2(:,2),'.','color',cmap(f,:));
            end
        catch er
            disp(['error on ',num2str(f)])
            disp(er.getReport);
            disp('stop here')
        end
    end
    toc
    disp('all points linked')
    
    % should check that there is extra data columns first?
    
    %%
    
    figure(11); clf;
    imagesc( squeeze(pntChain(:,1,:)) - squeeze(nanmean(pntChain(:,1,:),3))  ) ; colorbar; clim([-4 4]);
    GetColorMap('RedWhiteBlueSat');
    
    %% Pull the rest of the data into the point chains
    pntArray = nan(size(pntChain,1),nFrames,6);  % x,y,z,t,h,b 
    % (maybe this should be a structure for each spot?)  
    % (I kinda like matrices though
    % data(m(:,2),frameNum,all_data_types) = dataInFrame_f(m(:,1),all_data_types) 
    for f=1:nFrames  % f = 5
        m = idx{f};
        if ~isempty(m)
            pntArray(m(:,2),f,1) = fits1(f).x(m(:,1)); % that's compact but not super legible;  
            pntArray(m(:,2),f,2) = fits1(f).y(m(:,1)); % the pattern is helping with the reading
            pntArray(m(:,2),f,5) = fits1(f).height(m(:,1)); % spot brightness  (will need this later for calculating z)
            pntArray(m(:,2),f,4) = fits1(f).frame*ones(length(m(:,1)),1); % time (by frame)
            pntArray(m(:,2),f,6) = fits1(f).background(m(:,1)); % spot brightness  (will need this later for calculating z)
        end
    end
    
    % a quick look at the pnt chain data 
    figure(11); clf; ca = [-2,2];
    subplot(1,3,1); imagesc(  squeeze(pntArray(:,:,1)) -nanmean(squeeze(pntArray(:,:,1)),2) );  colorbar;   clim(ca); title('chain1 x')
    subplot(1,3,2); imagesc(  squeeze(pntArray(:,:,2)) -nanmean(squeeze(pntArray(:,:,2)),2) ); colorbar;   clim(ca);  title('chain1 y')
    subplot(1,3,3); imagesc(  squeeze(pntArray(:,:,5))  ); colorbar;  clim([0,600]);  title('chain1 brightness')
    xlabel('frame'); ylabel('spot ID'); 
    %%
    if ~pars.multiFit
        % Gaussian fit the brightness and read that in as z. 
        ftype = fittype('exp(-((x-mu_x)/(2*s)).^2)*a +b',...
                         'coeff', {'a','mu_x','b','s'},'ind',{'x'});
        [nPts,Zs,~] = size(pntArray);
    % converted to table output
    %   1 entry per spot
    fitsO = array2table(nan(nPts,12));
      fitsO.Properties.VariableNames={'x','y','xm','ym','h','z','zw','za','rsquare','zci_lower','zci_upper','ciratio'};
        for n=1:nPts  % n=77 
            if pars.verbose
                if rem(n,100)==0
                    disp(['fitting ',num2str(n/nPts*100,3),'% complete'])
                end
            end
            b = pntArray(n,:,5)' + pntArray(n,:,6)';
            b(isnan(b)) =min(pntArray(n,:,6));
            if isempty(pars.z_per_frame)
                z = (1:Zs)';  % I think we just want this to be offset value for the frame; 
            else
                z = pars.z_per_frame;
                if length(b) > length(z) % This should be unnecessary but by some evil some movies have more frames than recorded offset values.    
                    b(length(z)+1:end) = [];
                end
            end
            if sum(isnan(b)) == length(b)
                continue
            end
            
        
             try
               %    Gaussian fit
               [h,i] = max(b);
               sp = [max(b),    z(i),   median(b), pars.sigma];
               lb = [0,         0,      0,         pars.sigma/pars.sigmaRange]; 
               ub = [4*max(b),  max(z), max(b),    pars.sigma*pars.sigmaRange];
        
                [fitZ,gof] = fit(z,b,ftype,...
                     'StartPoint',sp,...
                     'Lower',lb,...
                     'Upper',ub,...
                     'TolFun',1e-6,'TolX',1e-6,'MaxIter',4000,'MaxFunEvals',600);
                 % figure(10); clf; plot(fitZ,z,b); pause(.1); xlabel('z'); ylabel('brightness');
        
                ci = confint(fitZ,.9);
                if (ci(2,2)-ci(1,2))/(max(z)-min(z)) > .5
                    poorSpot = ['detected poor fit on spot = ',num2str(n)];
                    if pars.verbose
                        disp(poorSpot);
                    end
                    if pars.showPoorFits
                        figure(10); clf; plot(fitZ,z,b); pause(.1); xlabel('z'); ylabel('brightness'); title(poorSpot); pause(.01);
                    end
                    % disp('place debug here')
                end
                if  rem(n,round(nPts/10)) == 0
                  figure(10); clf; plot(fitZ,z,b); pause(.1); xlabel('z'); ylabel('brightness'); title(['example fit on spot = ',num2str(n)]);
                  % pause;
                end
                 fitsO.x(n) = nanmean(pntArray(n,:,1));
                 fitsO.y(n) = nanmean(pntArray(n,:,2));
                 fitsO.xm(n) = pntArray(n,i,1);
                 fitsO.ym(n) = pntArray(n,i,2);
                 fitsO.h(n) = h;
                 fitsO.z(n) = fitZ.mu_x;
                 fitsO.zw(n) = fitZ.s;
                 fitsO.za(n) = fitZ.a;
                 fitsO.rsquare(n) = gof.rsquare;
                 fitsO.zci_lower(n) = ci(1,2);
                 fitsO.zci_upper(n) = ci(2,2);
                 fitsO.ciratio(n) = (ci(2,2)-ci(1,2))/(max(z)-min(z));
            catch er
                disp(er.getReport);
                disp('place debug here')
            end   
        end
    
    else
        %% Gaussian multifit
        % Gaussian fit the brightness and read that in as z. 
        ftype = fittype('exp(-((x-mu_x)/(2*s)).^2)*a +b',...
                         'coeff', {'a','mu_x','b','s'},'ind',{'x'});
        nn = 0; % counter for spots (may end up exceeding nPts);
        [nPts,Zs,~] = size(pntArray);
        fitsO = array2table(nan(nPts,12));
        fitsO.Properties.VariableNames={'x','y','xm','ym','h','z','zw','za','rsquare','zci_lower','zci_upper','ciratio'};
        for n=1:nPts  % n=77
            % ---- A little progress reporter
            if pars.verbose
                if rem(n,100)==0
                    disp(['fitting ',num2str(n/nPts*100,3),'% complete'])
                end
            end
            % ---- read the brightness and z-pos data out of the array
            b_dat = pntArray(n,:,5)' + pntArray(n,:,6)'; % brightness 
            b_dat(isnan(b_dat)) =min(pntArray(n,:,6));
            if isempty(pars.z_per_frame)
                z_pos = (1:Zs)';  % I think we just want this to be offset value for the frame; 
            else
                z_pos = pars.z_per_frame;
                if length(b_dat) > length(z_pos) % This should be unnecessary but by some evil some movies have more frames than recorded offset values.    
                    b_dat(length(z_pos)+1:end) = [];
                end
            end
            if sum(isnan(b_dat)) == length(b_dat)
                continue
            end
            % --- now we have the data "z_post" for z-pos, "b_dat" for brightness
           
            % === Multi fit processing
            try
                % ---- some key parameters are translated to shorthand
                N = length(b_dat);
                bkd = min(b_dat);
                res = pars.multiFit_res; 
                % ---- first we find multiple peaks along the trace:
                bf = b_dat'; % an editable copy of the data
                bf = bf - bkd;
                bf(bf<pars.multiFit_thetaAbs)  = 0;  % dim pixels to zero  (might not want to do this before resizeing) 
                bf = imresize(bf,[1,N/res]); % 
                % zds= imresize(z_pos',[1,N/res]);  % for troublshooting
                % figure(2); clf; plot(z_pos,b_dat,'.-'); hold on;  % for troublshooting
                % plot(zds,bf,'.-');                 % for troublshooting
                bf(bf<pars.multiFit_thetaMax*max(bf))  = 0;  
              %   plot(zds,bf,'.-'); legend('original','filt1','filt2');       % for troublshooting
                pk = regionprops(bf>0,'centroid'); % slower but better than peaks = find(zf > 0)*res; 
                pk = cat(1,pk.Centroid);
                if ~isempty(pk)
                peaks = pk(:,1)'*res;
                xs = round([peaks-pars.multiFit_zWindow*res; peaks+pars.multiFit_zWindow*res]); % indices
                xs(xs<1) = 1; xs(xs>length(z_pos)) = length(z_pos);       
                nZs = length(peaks);
                else
                    nZs = 0;
                end
                if  rem(n,round(nPts/10)) == 0
                    figure(10); clf;
                end
            catch er
                warning(er.getReport);
                disp(er.message);
                disp('place debug stop here')
            end
            try
                for k=1:nZs
                    nn=n+1;
                    z = z_pos(xs(1,k):xs(2,k));
                    b = b_dat(xs(1,k):xs(2,k));
                    
                    %    Gaussian fit
                    [h,i] = max(b);
                    sp = [max(b),    z(i),   median(b), pars.sigma];
                    lb = [0,         0,      0,         pars.sigma/pars.sigmaRange]; 
                    ub = [4*max(b),  max(z), max(b),    pars.sigma*pars.sigmaRange];
                    
                    [fitZ,gof] = fit(z,b,ftype,...
                         'StartPoint',sp,...
                         'Lower',lb,...
                         'Upper',ub,...
                         'TolFun',1e-6,'TolX',1e-6,'MaxIter',4000,'MaxFunEvals',600);
                     % figure(10); clf; plot(fitZ,z,b); pause(.1); xlabel('z'); ylabel('brightness');
                    
                     try
                    ci = confint(fitZ);
                    if (ci(2,2)-ci(1,2))/(max(z)-min(z)) > .5
                        poorSpot = ['detected poor fit on spot = ',num2str(n)];
                        if pars.verbose
                            disp(poorSpot);
                        end
                        % figure(10); clf; plot(fitZ,z,b); pause(.1); xlabel('z'); ylabel('brightness'); title(poorSpot);
                        % disp('place debug here')
                    end
                     catch er  % can't compute if too little data
                         ci=nan;
                     end
                              
                    if  rem(n,round(nPts/10)) == 0
                        zz = linspace(z(1),z(end),100);
                        yy = feval(fitZ,zz);
                        y0 = 1.2*ones(length(yy),1);
                        figure(10); 
                        plot(zz,yy,'.-','color',cmap(k,:)); hold on;
                        plot(zz,y0,'-','color',cmap(k,:)); hold on;
                        % figure(11); clf;  plot(zz,yy,'.-','color',cmap(k,:));
                    end
    
                     nn = nn+1;
                     fitsO.x(nn) = nanmean(pntArray(n,:,1));
                     fitsO.y(nn) = nanmean(pntArray(n,:,2));
                     fitsO.xm(nn) = pntArray(n,i,1);
                     fitsO.ym(nn) = pntArray(n,i,2);
                     fitsO.h(nn) = h;
                     fitsO.z(nn) = fitZ.mu_x;
                     fitsO.zw(nn) = fitZ.s;
                     fitsO.za(nn) = fitZ.a;
                     fitsO.rsquare(nn) = gof.rsquare;
                     fitsO.zci_lower(nn) = ci(1,2);
                     fitsO.zci_upper(nn) = ci(2,2);
                     fitsO.ciratio(nn) = (ci(2,2)-ci(1,2))/(max(z)-min(z));
                     % fitsO(nn).x = nanmean(pntArray(n,:,1));
                     % fitsO(nn).y = nanmean(pntArray(n,:,2));
                     % fitsO(nn).xm = pntArray(n,i,1);
                     % fitsO(nn).ym = pntArray(n,i,2);
                     % fitsO(nn).h = h;
                     % fitsO(nn).z = fitZ.mu_x;
                     % fitsO(nn).zw = fitZ.s;
                     % fitsO(nn).za = fitZ.a;
                     % fitsO(nn).rsquare = gof.rsquare;
                     % fitsO(nn).zci = ci(:,2); % ci just on z-pos
                     % fitsO(nn).ciratio = (ci(2,2)-ci(1,2))/(max(z)-min(z));
                end
                if  rem(n,round(nPts/10)) == 0
                      figure(10); plot(z_pos,b_dat,'k.-'); xlabel('z'); 
                      ylabel('brightness'); title(['example fit on spot = ',num2str(n)]);
                      pause(.1); 
                end
    
            catch er
                disp(er.getReport);
                disp('place debug here')
            end   
        end
    
    end
    
    
    if strcmp(pars.output,'struct') % for backwards compatibility
        fitsO = table2struct(fitsO);
    end

end