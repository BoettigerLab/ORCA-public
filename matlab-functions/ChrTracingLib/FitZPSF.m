function fitsO = FitZPSF(currTable,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'Zs','positive',[]};  % width
defaults(end+1,:) = {'sigma','positive',5};  % width
defaults(end+1,:) = {'sigmaRange','positive',4}; % fold change allowed for sigma.  Allowed widths will be = sigma/sigmaRange:sigma*sigmaRange
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'z_per_frame','array',[]}; % 
defaults(end+1,:) = {'figHandle','freeType',[]}; % 
defaults(end+1,:) = {'multiFit','boolean',false}; % 
defaults(end+1,:) = {'multiFit_res','positive',2}; % resolution (in pixels) at which to attempt distinguishing
defaults(end+1,:) = {'multiFit_thetaAbs','fraction',200}; %  threshold in absolute brightness units  
defaults(end+1,:) = {'multiFit_thetaMax','positive',.75}; % threshold relative to max brightness along the z-pixel stack
defaults(end+1,:) = {'multiFit_zWindow','positive',1.4}; % threshold in brightness units   
defaults(end+1,:) = {'showPoorFits','boolean',false}; %good for troubleshooting
defaults(end+1,:) = {'showEveryXFits','nonnegative',100}; %good for troubleshooting

pars = ParseVariableArguments(varargin,defaults,mfilename);


nPts = length(unique(currTable.traceID));

% Need to know how many z-steps 
%    if not given just assume range
if isempty(pars.Zs)
    Zs = max(currTable.frame) - min(currTable.frame) + 1;
else
    Zs = pars.Zs;
end


if ~pars.multiFit
    % Gaussian fit the brightness and read that in as z. 
    ftype = fittype('exp(-((x-mu_x)/(2*s)).^2)*a +b',...
                     'coeff', {'a','mu_x','b','s'},'ind',{'x'});

% converted to table output
%   1 entry per spot
fitsO = array2table(nan(nPts,13));
  fitsO.Properties.VariableNames={'x','y','z','obsPeak','zSigma','zPeak','background','rsquare','zci_lower','zci_upper','frames','nFrames','traceID'};

    for n=1:nPts  % n=77 
        traceTable = currTable(currTable.traceID==n,:);
        if isempty(traceTable)
            continue
        end
        if pars.verbose
            if rem(n,100)==0
                disp(['fitting ',num2str(n/nPts*100,3),'% complete'])
            end
        end
        try
        b = nan(Zs,1);
        fr = rem(traceTable.frame,Zs); fr(fr==0)=Zs; % convert frame position to frame height in cycle 
        b(fr) = traceTable.background + traceTable.height; % spot brightness
        weight = traceTable.height./sum(traceTable.height);
        xx = weight'*traceTable.x; % weighted sum, x position
        yy = weight'*traceTable.y; % weighted sum, x position
        bkd =  min(traceTable.background);
        b(isnan(b)) = bkd;
        catch er
            disp(er.getReport)
            disp('debug here')
        end
        
         z = (1:Zs)';  % Maybe we think we just want this to be offset value for the frame, provided we correct for the missing frames   
           
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
            if  rem(n,round(nPts/pars.showEveryXFits)) == 0 
              figure(10); clf; plot(fitZ,z,b); pause(.1); xlabel('z'); ylabel('brightness'); title(['example fit on spot = ',num2str(n)]);
              % pause;
            end
             fitsO.x(n) = xx;
             fitsO.y(n) = yy;
             fitsO.z(n) = fitZ.mu_x;
             fitsO.obsPeak(n) = h;
             fitsO.zSigma(n) = fitZ.s;
             fitsO.zPeak(n) = fitZ.a;
             fitsO.background(n) = bkd;
             fitsO.rsquare(n) = gof.rsquare;
             fitsO.zci_lower(n) = ci(1,2);
             fitsO.zci_upper(n) = ci(2,2);
             fitsO.ciratio(n) = (ci(2,2)-ci(1,2))/(max(z)-min(z));
             fitsO.frames(n) = mean(traceTable.frame);
             fitsO.nFrames(n) = length(traceTable.frame);
             fitsO.traceID(n) = n;
             
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