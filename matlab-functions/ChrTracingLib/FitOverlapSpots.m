function  [fitResults,ci95] = FitOverlapSpots(im,nEmits,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'method',{'fit','lsqn'},'lsqn'};
defaults(end+1,:) = {'dimensions',{'2D','3D'},'2D'};
defaults(end+1,:) = {'troubleshoot','boolean',false};
defaults(end+1,:) = {'showplot','boolean',true};
%  ROI selection (not active yet, should be passed to sep function?) 
% defaults(end+1,:) = {'peakBlur','nonnegative',0};
% defaults(end+1,:) = {'keepBrightest','integer',inf};  % coud use more general  
% defaults(end+1,:) = {'trimBorder','nonnegative',0};
% defaults(end+1,:) = {'cropBoxWidth','nonnegative',8};
% fitter parameters
defaults(end+1,:) = {'FiniteDifferenceType',{'central','forward'},'central'};
defaults(end+1,:) = {'MaxFunctionEvaluations','positive',5E4};
defaults(end+1,:) = {'OptimalityTolerance','positive',1E-9};
% units
defaults(end+1,:) = {'xyUnitConvert','positive',1}; % Not active yet
% PSF initialization and bounds
defaults(end+1,:) = {'initXY','freeType',[]}; % empty for auto compute
defaults(end+1,:) = {'initPeak','freeType',[]}; % empty for auto compute
defaults(end+1,:) = {'initBkd','freeType',[]};% empty for auto compute
defaults(end+1,:) = {'initSigma','nonnegative',1.4};
defaults(end+1,:) = {'minPeakHeight','nonnegative',1};

defaults(end+1,:) = {'minSigma','nonnegative',.1};
defaults(end+1,:) = {'maxSigma','nonnegative',2.5};
defaults(end+1,:) = {'minHeight','nonnegative',0};
defaults(end+1,:) = {'maxHeight','nonnegative',inf};
defaults(end+1,:) = {'maxShift','nonnegative',inf};
defaults(end+1,:) = {'cameraBackground','nonnegative',0};
pars = ParseVariableArguments(varargin,defaults,mfilename);


% nEmits = length(xy);

data_2d = double(im - pars.cameraBackground);
[rows,cols] = size(data_2d);
[X,Y] = meshgrid(1:cols, 1:rows); %#ok<*ASGLU>  % USED by eval expression

%===== if no start points are passed ======
if isempty(pars.initXY)
    bw = imregionalmax(data_2d); %  figure(2); clf; imagesc(bw);
    bw(data_2d<pars.minPeakHeight) = 0;
    localPeaks = find(bw(:));
    [mu_y0,mu_x0] = ind2sub(size(bw), localPeaks);
else
    mu_x0 = pars.initXY(:,1);
    mu_y0 = pars.initXY(:,2);
end
    
% ==== default PSF properties (should be parameters) ====
if isempty(pars.initPeak)
    % a0 = max(data_2d(:));
    a0 = zeros(length(mu_y0),1);
    for n=1:length(a0)
        a0(n) = data_2d(round(mu_y0(n)),round(mu_x0(n)));
    end
else
    a0 = pars.initPeak; 
end
if isempty(pars.initBkd)
    b0 = min(data_2d(:));
else
    b0 = pars.initBkd;
end
sigma_x0 = pars.initSigma;
sigma_y0 = pars.initSigma; % enforce symmetric for now
minSigma = pars.minSigma;
maxSigma = pars.maxSigma;

    %--------------------------------------------------------------------%  
    % least squares fitting to a sum of N 2D Gaussians
    %    see FitPsf3D.m for generalization to 3D using lsqn fitting
%     par0 = [0,sigma_x0,0,sigma_y0,a0,b0]; % initial fits, all N gaussians in same place 
%     par0 = repmat(par0,1,nEmits);
    par0 = [sigma_x0,b0,zeros(1,nEmits*3)];
    % if number of emitters doesn't match number of starting peaks
    %  (should only happen if starting peaks are automatically determined)
    %  If too few, just replicate. If too many take the most separated?
    nlmax = nEmits-length(mu_x0);
    if nlmax > 0
       mu_y0 = [mu_y0; repmat(mu_y0(1),nlmax,1)];
       mu_x0 = [mu_x0; repmat(mu_x0(1),nlmax,1)];
       a0    = [a0; repmat(a0(1),nlmax,1)];
    elseif nlmax <0  
        mu_x0 = mu_x0(1:nEmits);
        mu_y0 = mu_y0(1:nEmits); 
        a0    = a0(1:nEmits);
    end
    for n=1:nEmits
        % par0(1+(n-1)*6) = mu_x0(n);
        % par0(3+(n-1)*6) = mu_y0(n);
        par0(2+1+(n-1)*3) = mu_x0(n);
        par0(2+2+(n-1)*3) = mu_y0(n);
        par0(2+3+(n-1)*3) = a0(n);
    end
      
    switch pars.method
        case 'lsqn'
            gaus = 'exp(-((X(:)-muX)/(2*sXY)).^2 -((Y(:)-muY)/(2*sXY)).^2 )*a +';
            gaus = regexprep(gaus,{'sXY'},{'p(1)'});
            % enforces symmetric XY right now
            %   should allow deformations to correct for comma aberration. 
            
            emitPars = cell(nEmits,1);
            emitExpr = cell(nEmits,1);
            for n=1:nEmits
                  emitPars{n} = {['p(',num2str(2+1+(n-1)*3),')'],['p(',num2str(2+2+(n-1)*3),')'],['p(',num2str(2+3+(n-1)*3),')']};
                  emitExpr{n} =  regexprep(gaus,{'muX','muY','a'},emitPars{n});
            end
            mGauss = cat(2,emitExpr{:});
            expr = ['minRes = @(p)', mGauss, '+p(3)  -data_2d(:);'];
            eval(expr);   

            % upper and lower bounds
            lb = [minSigma,0,zeros(1,nEmits*3)];
            lb(2+1:3:end)= mu_x0 - pars.maxShift;
            lb(2+2:3:end)= mu_y0 - pars.maxShift;
            lb(2+3:3:end) = pars.minHeight;
            ub = [maxSigma,inf,inf*ones(1,nEmits*3)];
            ub(2+1:3:end)= mu_x0 + pars.maxShift;
            ub(2+2:3:end)= mu_y0 + pars.maxShift;
            ub(2+3:3:end) = pars.maxHeight;
%             lb = [-inf,minSigma,-inf,minSigma,0,0]; % lower bound 
%             lb = repmat(lb,1,nEmits);
%             lb(1:6:6*(nEmits-1)+1) = mu_x0 - pars.maxShift;
%             lb(3:6:6*(nEmits-1)+3) = mu_y0 - pars.maxShift;
%             ub = [inf,maxSigma,inf,maxSigma,inf,inf]; % upper bound
%             ub = repmat(ub,1,nEmits);
%             ub(1:6:6*(nEmits-1)+1) = mu_x0 + pars.maxShift;
%             ub(3:6:6*(nEmits-1)+3) = mu_y0 + pars.maxShift;

            options = optimoptions('lsqnonlin',...
                                   'MaxFunctionEvaluations',pars.MaxFunctionEvaluations,...
                                   'OptimalityTolerance',pars.OptimalityTolerance,...
                                   'Display','off');

            [par,resnorm,residual,~,~,~,jacobian] = lsqnonlin(minRes, par0,lb,ub,options);   
            ci95 = nlparci(par,residual,'jacobian',jacobian,'alpha',1-.95);
   
 %%   
        case 'fit'
            gaus = 'exp(-((x-muX)/(2*sX)).^2 -((v-muY)/(2*sX)).^2 )*a0+';
            emitPars = cell(nEmits,1);
            emitExpr = cell(nEmits,1);
            for n=1:nEmits
                emitPars{n} = {['mx',num2str(n)],['my',num2str(n)],['a0',num2str(n)]};
                emitExpr{n} =  regexprep(gaus,{'muX','muY','a0'},emitPars{n});
                % emitPars{n} = {['mx',num2str(n)],['sx',num2str(n)],['my',num2str(n)],['sy',num2str(n)],['a',num2str(n)],['b',num2str(n)]};
                % emitExpr{n} =  regexprep(gaus,{'muX','sX','muY','sY','a','b'},emitPars{n});
            end
            mGauss = cat(2,'b+',emitExpr{:});
            mPars = cat(2,'sx','b',emitPars{:});
            mGauss = mGauss(1:end-1);
            ftype = fittype(mGauss,'coeff', mPars,'ind',{'x','v'});
            %  ftype = fittype('exp(-((x-mu_x)/(2*sigma_x)).^2-((v-mu_y)/(2*sigma_y)).^2 )*a/(2*pi*sigma_x*sigma_y)+b',...
            %                         'coeff', {'a','mu_x','sigma_x','mu_y','sigma_y','b'},'ind',{'x','v'}); 
            %    % dataTest = a0*exp(-((X-mu_x0-1.3)/(2*sigma_x0+3)).^2-((Y-mu_y0+.25)/(2*sigma_y0+.5)).^2 )+b0;
            fit_2d = fit([Y(:),X(:)],double(data_2d(:)),ftype,'StartPoint',par0);
            ci95 = confint(fit_2d);
            par = coeffvalues(fit_2d);      
    end
    fitResults.x = par(2+1:3:end);
    fitResults.y = par(2+2:3:end);
    fitResults.a = par(2+3:3:end);
    fitResults.s = par(1);
    fitResults.b = par(2);
    fitResults.ci95 = ci95;

    imVal = minRes(par) + data_2d(:);
    imVal = reshape(imVal,rows,cols);

    fitResults.imEst = imVal;

    if pars.showplot
     figure(10); clf;
     subplot(1,2,1); imagesc(data_2d); colorbar; hold on; plot(fitResults.x,fitResults.y,'ro'); title('data xy');
     subplot(1,2,2); imagesc(imVal); colorbar;  hold on; plot(fitResults.x,fitResults.y,'ro'); title('fit xy');
     colormap(gray);
    end
    