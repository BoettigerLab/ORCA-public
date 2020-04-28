function  [fitResults,ci95] = FitOverlapSpots3D(im,nEmits,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'method',{'fit','lsqn'},'lsqn'};
defaults(end+1,:) = {'dimensions',{'2D','3D'},'3D'}; % no active yet, could generalize later
defaults(end+1,:) = {'troubleshoot','boolean',false};
defaults(end+1,:) = {'showplot','boolean',true};
%  ROI selection - sent to FindPeaks3D
defaults(end+1,:) = {'peakBlur', 'nonnegative', .5}; % gaussian smoothing before initial find max 
defaults(end+1,:) = {'maxFitWidth', 'positive', inf}; % crop
defaults(end+1,:) = {'maxFitZdepth', 'positive', inf}; % crop
defaults(end+1,:) = {'keepBrightest','integer',inf};
defaults(end+1,:) = {'relativeHeight','fraction',0};
defaults(end+1,:) = {'minSep','positive',0}; % minimum separation in pixels between initial maxima
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'veryVerbose','boolean',false};
defaults(end+1,:) = {'troubleshoot','boolean',true};
% fitter parameters
defaults(end+1,:) = {'FiniteDifferenceType',{'central','forward'},'central'};
defaults(end+1,:) = {'MaxFunctionEvaluations','positive',5E4};
defaults(end+1,:) = {'OptimalityTolerance','positive',1E-9};
% units
defaults(end+1,:) = {'xyUnitConvert','positive',1}; % Not active yet
defaults(end+1,:) = {'zUnitConvert','positive',1}; % Not active yet
% PSF initialization and bounds
defaults(end+1,:) = {'initXYZ','freeType',[]}; % empty for auto compute
defaults(end+1,:) = {'initPeak','freeType',[]}; % empty for auto compute
defaults(end+1,:) = {'initBkd','freeType',[]};% empty for auto compute
defaults(end+1,:) = {'initSigma','nonnegative',1.4};
defaults(end+1,:) = {'initSigmaZ','nonnegative',2};
defaults(end+1,:) = {'minPeakHeight','nonnegative',1};
% bounds for fitting
defaults(end+1,:) = {'minSigma','nonnegative',.1};
defaults(end+1,:) = {'maxSigma','nonnegative',2.5};
defaults(end+1,:) = {'minSigmaZ','nonnegative',1};
defaults(end+1,:) = {'maxSigmaZ','nonnegative',10};
defaults(end+1,:) = {'minHeight','nonnegative',0};
defaults(end+1,:) = {'maxHeight','nonnegative',inf};
defaults(end+1,:) = {'maxShift','nonnegative',inf};
defaults(end+1,:) = {'minBkd','nonnegative',0};
defaults(end+1,:) = {'maxBkd','nonnegative',inf};
defaults(end+1,:) = {'cameraBackground','nonnegative',0};
pars = ParseVariableArguments(varargin,defaults,mfilename);

%% main function
% 
data_3d = double(im - pars.cameraBackground);
[rows,cols,stcks] = size(data_3d);
[Y,X,Z] = meshgrid(1:cols, 1:rows, 1:stcks); %#ok<ASGLU>  these are used.  

% Find peaks
if isempty(pars.initXYZ) 
   dTable = FindPeaks3D(im,'parameters',pars,'troubleshoot',false,'keepBrightest',nEmits);
else
    nPts = size(pars.initXYZ,1);
    x = zeros(nPts,1);
    y = zeros(nPts,1);
    z = zeros(nPts,1);
    h = zeros(nPts,1);
    for n=1:nPts
        x(n) = pars.initXYZ(:,1);
        y(n) = pars.initXYZ(:,2);
        z(n) = pars.initXYZ(:,3);
        h(n) = data_3d(x,y,z); 
    end
    dTable = table(x,y,z,h);
end
nEmitFound = height(dTable); 
if isempty(pars.initBkd)
    b0 = min(data_3d(:));
else
    b0 = pars.initBkd;
end
    %--------------------------------------------------------------------%  
    % least squares fitting to a sum of N 3D Gaussians
    par0 = [pars.initSigma,pars.initSigmaZ,b0,zeros(1,nEmits*4)];
    nlmax = nEmits-nEmitFound;
    if nlmax > 0
        while nEmitFound < nEmits
            dTable{end+1,:} = dTable{end,:} + rand( size(dTable(end,:)) ); %#ok<AGROW>
            nEmitFound = height(dTable); 
        end    
    elseif nlmax <0  
        dTable = dTable(1:nEmits,:);
    end
    for n=1:nEmits
        par0(3+1+(n-1)*4) = dTable.x(n);
        par0(3+2+(n-1)*4) = dTable.y(n);
        par0(3+3+(n-1)*4) = dTable.z(n);
        par0(3+4+(n-1)*4) = dTable.h(n);
    end
      
    switch pars.method
        case 'lsqn'
            gaus = 'exp(-((Y(:)-muX)/(2*sXY)).^2 -((X(:)-muY)/(2*sXY)).^2  -((Z(:)-muZ)/(2*sZ)).^2 )*a +';  % following X,Y swap in FitPsf3D
            gaus = regexprep(gaus,{'sXY','sZ'},{'p(1)','p(2)'});
            % parameters: sXY, sZ, b0, muX_1, muY_1, muZ_1, a_1, muX_2,...     
            emitPars = cell(nEmits,1);
            emitExpr = cell(nEmits,1);
            for n=1:nEmits
                  emitPars{n} = {['p(',num2str(3+1+(n-1)*4),')'],['p(',num2str(3+2+(n-1)*4),')'],...
                                 ['p(',num2str(3+3+(n-1)*4),')'],['p(',num2str(3+4+(n-1)*4),')']};
                  emitExpr{n} =  regexprep(gaus,{'muX','muY','muZ','a'},emitPars{n});
            end
            mGauss = cat(2,emitExpr{:});
            expr = ['minRes = @(p)', mGauss, ' p(3)  - data_3d(:);'];
            eval(expr);   
            % upper and lower bounds
            lb = [pars.minSigma,pars.minSigmaZ,pars.minBkd,zeros(1,nEmits*4)];
            ub = [pars.maxSigma,pars.maxSigmaZ,pars.maxBkd,inf*ones(1,nEmits*4)];
            for n=1:nEmits
                lb(3+1+(n-1)*4) = dTable.x(n)- pars.maxShift;
                lb(3+2+(n-1)*4) = dTable.y(n)- pars.maxShift;
                lb(3+3+(n-1)*4) = dTable.z(n)- pars.maxShift;
                lb(3+4+(n-1)*4) = pars.minHeight;
                ub(3+1+(n-1)*4) = dTable.x(n)+ pars.maxShift;
                ub(3+2+(n-1)*4) = dTable.y(n)+ pars.maxShift;
                ub(3+3+(n-1)*4) = dTable.z(n)+ pars.maxShift;
                ub(3+4+(n-1)*4) = pars.maxHeight;
            end
            options = optimoptions('lsqnonlin',...
                                   'MaxFunctionEvaluations',pars.MaxFunctionEvaluations,...
                                   'OptimalityTolerance',pars.OptimalityTolerance,...
                                   'Display','off');
            [par,resnorm,residual,~,~,~,jacobian] = lsqnonlin(minRes, par0,lb,ub,options);   
            ci95 = nlparci(par,residual,'jacobian',jacobian,'alpha',1-.95); 
 %%   
        case 'fit'
            error('fxn "fit" does not work in 3D in matlab as of R2018');
            gaus = 'exp(-((x-muX)/(2*sXY)).^2 -((v-muY)/(2*sXY)).^2  -((z-muZ)/(2*sZ)).^2 )*a+'; %#ok<*UNRCH>
            emitPars = cell(nEmits,1);
            emitExpr = cell(nEmits,1);
            for n=1:nEmits
                emitPars{n} = {['mx',num2str(n)],['my',num2str(n)],['mz',num2str(n)],['a',num2str(n)]};
                emitExpr{n} =  regexprep(gaus,{'muX','muY','muZ','a'},emitPars{n});
            end
            mGauss = cat(2,'b+',emitExpr{:});
            mPars = cat(2,'sXY','sZ','b',emitPars{:});
            mGauss = mGauss(1:end-1);
            ftype = fittype(mGauss,'coeff', mPars,'ind',{'x','v','z'});
            fit_2d = fit([Y(:),X(:),Z(:)],double(data_3d(:)),ftype,'StartPoint',par0);
            ci95 = confint(fit_2d);
            par = coeffvalues(fit_2d);      
    end
    %% export results in a more readable format
    fitResults.x = par(3+1:4:end);
    fitResults.y = par(3+2:4:end);
    fitResults.z = par(3+3:4:end);
    fitResults.a = par(3+4:4:end);
    fitResults.sXY = par(1);
    fitResults.sZ = par(2);
    fitResults.b = par(3);
    fitResults.dTable = dTable;
    % fitResults.ci95 = ci95;
    imVal = minRes(par)+ data_3d(:);
    imVal = reshape(imVal,rows,cols,stcks);
    fitResults.imEst = imVal;
    if pars.showplot
        figure(10); clf;
        subplot(2,2,1); imagesc(max(data_3d,[],3)); colorbar; hold on; plot(fitResults.x,fitResults.y,'ro');%  title('data xy');
        subplot(2,2,2); imagesc(max(imVal  ,[],3)); colorbar; hold on; plot(fitResults.x,fitResults.y,'ro'); title('fit xy');
        subplot(2,2,3); imagesc(max(permute(data_3d,[3,2,1]),[],3));  colorbar;  hold on; plot(fitResults.x,fitResults.z,'ro');%  title('data xz');
        subplot(2,2,4); imagesc(max(permute(imVal,[3,2,1]),[],3));   colorbar;  hold on; plot(fitResults.x,fitResults.z,'ro');  title('fit xz');
        colormap(gray);
    end
    