function  [fitResults,fitError] = FitOverlapSpots3D(im,nEmits,varargin)
% Finds 'nEmits' emitters in a single image, allowing for overlapping.
%   this functions solves the full least squares minimization for N total
%   emitters lying in the same 3D region.  Proper bounds on the PSFs,
%   especially the relative brightness of these hidden emitters is needed.
% Critically, the number of emitters is assumed to be known. For variable
%   numbers of emitters, you can relax the bounds (specifically the height
%   variation restriction 'brightVar'). 
% 
% to do ??
% 
% An explicit, stochastic noise model.
% we don't want to use the PSF centers/width to capture features of the
% image that are due to the pixel noise. 
% I believe this is the Weiner filter
% 
% Filter for 'spot like' properties: if one of the spots does not meet the
% required filters, exclude it and use instead the brightest remaining spot
% Maybe this is best done in downstream processing by a different function?
% 
% For future: 
% I suppose a better way is to compute both a noise-only fit, and 1-fit,
% 2-fit up to N-fit and decide which is best, with forcing N fits as an
% option. 

defaults = cell(0,3);
defaults(end+1,:) = {'method',{'fit','lsqn'},'lsqn'};
defaults(end+1,:) = {'dimensions',{'2D','3D'},'3D'}; % no active yet, could generalize later
defaults(end+1,:) = {'troubleshoot','boolean',false};
defaults(end+1,:) = {'showplot','boolean',false};
defaults(end+1,:) = {'figShowFits','freeType',11};
defaults(end+1,:) = {'reconstruct','boolean',false}; % reconstruct a noise-free image from the fit values to compare to the original data. 
% PSF initialization and bounds and ROI detection
defaults(end+1,:) = {'minDif','nonnegative',.5}; % extra peaks must be within this factor of the bright peak
defaults(end+1,:) = {'initXYZ','freeType',[]}; % empty for auto compute
defaults(end+1,:) = {'initPeak','freeType',[]}; % empty for auto compute
defaults(end+1,:) = {'initBkd','freeType',[]};% empty for auto compute
defaults(end+1,:) = {'initSigma','nonnegative',1.4};
defaults(end+1,:) = {'initSigmaZ','nonnegative',3};
defaults(end+1,:) = {'brightVar','nonnegative',1/5};  % 1/5 is quite strict for equally bright. 10 or inf is very relaxed  
% fitter parameters
defaults(end+1,:) = {'FiniteDifferenceType',{'central','forward'},'central'};
defaults(end+1,:) = {'MaxFunctionEvaluations','positive',5E4};
defaults(end+1,:) = {'OptimalityTolerance','positive',1E-10};
% units
%    ??  maybe these should always be handled outside of the fitter? in
%    principle the PSF is related to the units too, and those values are in
%    pixels. 
defaults(end+1,:) = {'xyUnitConvert','positive',1}; % Not active yet  
defaults(end+1,:) = {'zUnitConvert','positive',1}; % Not active yet
% bounds for fitting
defaults(end+1,:) = {'minHeight','nonnegative',0}; % could be a vector in descending order if different heights are expected  
defaults(end+1,:) = {'maxHeight','nonnegative',inf}; % could be a vector in descending order if different heights are expected
defaults(end+1,:) = {'minSigma','nonnegative',.1}; % scalar. default is generous, multifitting will be more accurate if properly constrained
defaults(end+1,:) = {'maxSigma','nonnegative',2.5}; % scalar  default is generous, multifitting will be more accurate if properly constrained
defaults(end+1,:) = {'minSigmaZ','nonnegative',.5}; % scalar default is generous, multifitting will be more accurate if properly constrained
defaults(end+1,:) = {'maxSigmaZ','nonnegative',10}; % scalar default is generous, multifitting will be more accurate if properly constrained
defaults(end+1,:) = {'maxShift','nonnegative',inf};
defaults(end+1,:) = {'minBkd','nonnegative',0};
defaults(end+1,:) = {'maxBkd','nonnegative',inf};
defaults(end+1,:) = {'maxTilt','nonnegative',0}; % a slanted plane background model gaus(X,Y) +px*X +py*Y +b0.  maxTilt sets the amplitude of px and py.
defaults(end+1,:) = {'cameraBackground','nonnegative',0};
pars = ParseVariableArguments(varargin,defaults,mfilename);

%% main function
% 
if pars.initSigmaZ < pars.minSigmaZ
    warning('invalid: pars.initSigmaZ < pars.minSigmaZ');
    pars.initSigmaZ = pars.minSigmaZ;
elseif pars.initSigmaZ > pars.maxSigmaZ
    warning('invalid: pars.initSigmaZ > pars.maxSigmaZ');
    pars.initSigmaZ = pars.maxSigmaZ;
end


data_3d = double(im - pars.cameraBackground);
% Find peaks
if isempty(pars.initXYZ) 
    % find regional maxima that are within 'minDif' from one-another
    %   could expand into its own function with additional optional pars
    %   keepBrightest, minSep, peak blur etc.  
    % 
    % data_3d = datSpot
    imagePeaks = data_3d;
    medV = median(imagePeaks(:));
    imagePeaks = imagePeaks - medV;
    maxV = max(imagePeaks(:));
    imagePeaks(imagePeaks<pars.minDif*maxV) = 0;
    bw = imregionalmax(imagePeaks);
    [y,x,z] = ind2sub(size(imagePeaks),find(bw));
    h = data_3d(bw);
    [~,i] = sort(h,'descend');
    dTable = table(x,y,z,h);
    dTable = dTable(i,:); % sort
    
    if height(dTable) < nEmits
        % Just take the brightest N points in the 3D image as staring foci
        [~,id] = sort(data_3d(:),'descend') ; % data_3d = im3D;
        x = zeros(nEmits,1);
        y = zeros(nEmits,1);
        z = zeros(nEmits,1);
        h = zeros(nEmits,1);
        for n=1:nEmits
            [y(n),x(n),z(n)] = ind2sub(size(data_3d),id(n));
            h(n) = data_3d(y(n),x(n),z(n));
        end   
        dTable = table(x,y,z,h);
    end
else
    nPts = size(pars.initXYZ,1);
    x = zeros(nPts,1);
    y = zeros(nPts,1);
    z = zeros(nPts,1);
    h = zeros(nPts,1);
    for n=1:nPts
        x(n) = pars.initXYZ(n,1);
        y(n) = pars.initXYZ(n,2);
        z(n) = pars.initXYZ(n,3);
        h(n) = data_3d(round(y(n)),round(x(n)),round(z(n)));  % rows by cols
    end
    dTable = table(x,y,z,h);
end
nEmitFound = height(dTable); 

% compress image size to improve fitting speed by trimming out extra pixels
%   thus, is a large image is passed, but all the spots are in the center,
%   the fitting will go faster as the data sent to the least squares
%   optimization is much smaller. But if the whole volume is scattered with
%   data the whold volume will be used.
[ymax,xmax,zmax] = size(data_3d);
x1 = max(1, floor(min(x) - 2*pars.maxSigma));
x2 = min(xmax, floor(max(x) + 2*pars.maxSigma));
y1 = max(1, floor(min(y) - 2*pars.maxSigma));
y2 = min(ymax, floor(max(y) + 2*pars.maxSigma));
z1 = max(1, floor(min(z) - 2*pars.maxSigmaZ));
z2 = min(zmax, floor(max(z) + 2*pars.maxSigmaZ));
data_3dcrop = data_3d(y1:y2,x1:x2,z1:z2);

if isempty(pars.initBkd)
    b0 = median(data_3dcrop(:));
else
    b0 = pars.initBkd;
end

% brightness estimates from data
% Here we get the max and min brightness from the data
%   we can't simply use the brightest pixel as the max, sometimes the
%   brightest pixel will be made of N stacked emitters, but sometimes the N
%   emitters are all separate.  
%   Here, we look at the separation of the two brightest emitters to first
%   decide if stacking has occurred and too what degree (based on the
%   expected PSF decay for this separation).  This tells us how much 
%   This should be corrected for the case where more than 2 emitters stack.
if height(dTable) > 1
    initDists = pdist(dTable{1:2,1:3}); % separation between the two brightest spots (could be extended)
    mult = exp(-(initDists./(2*pars.initSigma)).^2) ; % expected gaussian decay for this separation. this = 1 if separation = 0, and the two spots are stacked
    aveBright = (mean(dTable.h(1:2))-b0)/(1+mult); % 
    aveBright = max([0,aveBright]);
else
    aveBright = dTable.h;
end

% how about a 2D fit of a single xy plane and a single xz plane, where the
% two planes are selected to pass through the brightest pixel? 
%   This could be faster tor single-fitters, it doesn't generalize to
%   multifitting though...

if pars.minHeight == 0  % height is Gaussian peak height (i.e. max minus background) 
    for n=1:nEmits
       pars.minHeight(n) = (1-pars.brightVar)*aveBright;
    end
end
if isinf(pars.maxHeight)
    for n=1:nEmits
       pars.maxHeight(n) = (1+pars.brightVar)*aveBright;
    end
end

% match length
%    especially when minHeight maxHeight are measured but assumed identical
while length(pars.minHeight)<nEmits
    pars.minHeight(end+1) = pars.minHeight(end);
end
while length(pars.maxHeight)<nEmits
    pars.maxHeight(end+1) = pars.maxHeight(end);
end

    %--------------------------------------------------------------------%  
    % least squares fitting to a sum of N 3D Gaussians
    
    [rows,cols,stcks] = size(data_3dcrop);
    [X,Y,Z] = meshgrid(1:cols, 1:rows, 1:stcks); %#ok<ASGLU> % updated match meshgrid.  
    
    
    % if fewer initial points are specified than the N centroids requested,
    % initialize new points with random shifts from last starting point.
    nlmax = nEmits-nEmitFound;
    if nlmax > 0
        while nEmitFound < nEmits
            dTable{end+1,:} = dTable{end,:} + .5*rand( size(dTable(end,:)) ); %#ok<AGROW>
            nEmitFound = height(dTable); 
        end    
    elseif nlmax <0  
        dTable = dTable(1:nEmits,:);
    end
    
    % shared pars  sXY, sZ, b, px, py,
    commonPars = [pars.initSigma,pars.initSigmaZ,b0,0,0];
    f = length(commonPars); 
    par0 = [commonPars,zeros(1,nEmits*4)];
    for n=1:nEmits
        par0(f+1+(n-1)*4) = dTable.x(n);
        par0(f+2+(n-1)*4) = dTable.y(n);
        par0(f+3+(n-1)*4) = dTable.z(n);
        par0(f+4+(n-1)*4) = max([0,dTable.h(n)-b0]);
    end
      
    fitSuccess = true;
    
    switch pars.method
        case 'lsqn'
            gaus = 'exp(-((X(:)-muX)/(2*sXY)).^2 -((Y(:)-muY)/(2*sXY)).^2  -((Z(:)-muZ)/(2*sZ)).^2 )*a +';  
            gaus = regexprep(gaus,{'sXY','sZ'},{'p(1)','p(2)'});
            % construct a string for the sum of N Gaussians...
            % parameters: sXY, sZ, b0, muX_1, muY_1, muZ_1, a_1, muX_2,...     
            emitPars = cell(nEmits,1);
            emitExpr = cell(nEmits,1);
            emitVars = cell(4,nEmits);
            for n=1:nEmits
                  emitPars{n} = {['p(',num2str(f+1+(n-1)*4),')'],['p(',num2str(f+2+(n-1)*4),')'],...
                                 ['p(',num2str(f+3+(n-1)*4),')'],['p(',num2str(f+4+(n-1)*4),')']};
                  emitExpr{n} =  regexprep(gaus,{'muX','muY','muZ','a'},emitPars{n});
                  emitVars{n} = {['muX',num2str(n)],['muY',num2str(n)],['muZ',num2str(n)],['a',num2str(n)]};
            end
            mGauss = cat(2,emitExpr{:});
            expr = ['minRes = @(p)', mGauss, ' p(3) +p(4)*X(:) +p(5)*Y(:)  - data_3dcrop(:);'];
            eval(expr);    % just builds the anonymous function to pass to lsqn
            % upper and lower bounds
            lb = [pars.minSigma,pars.minSigmaZ,pars.minBkd,-pars.maxTilt,-pars.maxTilt,zeros(1,nEmits*4)];
            ub = [pars.maxSigma,pars.maxSigmaZ,pars.maxBkd,pars.maxTilt,pars.maxTilt,inf*ones(1,nEmits*4)];
            for n=1:nEmits
                lb(f+1+(n-1)*4) = min([max([dTable.x(n)- pars.maxShift,1]),cols]);
                lb(f+2+(n-1)*4) = min([max([dTable.y(n)- pars.maxShift,1]),rows]);
                lb(f+3+(n-1)*4) = min([max([dTable.z(n)- pars.maxShift,1]),stcks]);
                lb(f+4+(n-1)*4) = pars.minHeight(n);
                ub(f+1+(n-1)*4) =  min([max([dTable.x(n)+ pars.maxShift,1]),cols]);
                ub(f+2+(n-1)*4) =  min([max([dTable.y(n)+ pars.maxShift,1]),rows]);
                ub(f+3+(n-1)*4) =  min([max([dTable.z(n)+ pars.maxShift,1]),stcks]);
                ub(f+4+(n-1)*4) = pars.maxHeight(n);
            end
            options = optimoptions('lsqnonlin',...
                                   'MaxFunctionEvaluations',pars.MaxFunctionEvaluations,...
                                   'OptimalityTolerance',pars.OptimalityTolerance,...
                                   'Display','off');
            if nargout < 2
                par = lsqnonlin(minRes,par0,lb,ub,options);   % Here we call the non-linear least squares fitting 
            else
                [par,~,residual,~,~,~,jacobian] = lsqnonlin(minRes,par0,lb,ub,options);   % This version also returns the error-stats non-linear least squares fitting
                try
                    ci95 = nlparci(par,residual,'jacobian',jacobian,'alpha',1-.95); 
                catch
                    ci95 = [];
                end
            end
            if sum(abs(par0-par))==0
                fitSuccess = false;
                warning('output = input');
            end
 %%   
        case 'fit'
            % Keep this for now
            % useful still for the collapse back to 2D version
            % Maybe future version of matlab will generlaize fit to 3D
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
            fit_2d = fit([X(:),Y(:),Z(:)],double(data_3dcrop(:)),ftype,'StartPoint',par0);
            ci95 = confint(fit_2d);
            par = coeffvalues(fit_2d);      
    end
    %% export results in a more readable format
    fitResults.x = par(f+1:4:end)+(x1-1);
    fitResults.y = par(f+2:4:end)+(y1-1);
    fitResults.z = par(f+3:4:end)+(z1-1);
    fitResults.a = par(f+4:4:end);
    fitResults.sXY = par(1);
    fitResults.sZ = par(2);
    fitResults.b = par(3);
    fitResults.px = par(4);
    fitResults.py = par(5);
    fitResults.fitQuality = fitSuccess;
    % these are for convience comparing
    fitResults.par = par;
    fitResults.par0 = par0;
    fitResults.parL = lb;
    fitResults.parU = ub;
    fitResults.dTable = dTable;
    fitResults.parNames = cat(2,{'sXY','sZ','b','px','py'},emitVars{:});
    if nargout > 1
        fitError.ci95 = ci95;
        fitError.residual = residual;
    end
    if pars.reconstruct
        imVal = minRes(par)+ double(data_3dcrop(:));  % imVal =  double(data_3d(:));
        imVal = reshape(imVal,rows,cols,stcks);
        fitResults.imEst = imVal;
    else
        fitResults.imEst = 0;
    end
    if pars.showplot && pars.figShowFits
        if pars.reconstruct
            figure(pars.figShowFits); clf;
            [xy,xz] = ProjectIm3D(data_3d); % use the original data
            [xyS,xzS] = ProjectIm3D(imVal);
            subplot(2,2,1); imagesc(xy);  colorbar; hold on; plot(fitResults.x,fitResults.y,'ro'); title('data xy');
            subplot(2,2,2); imagesc(xyS); colorbar; hold on; plot(par(f+1:4:end),par(f+2:4:end),'ro'); title('fit xy');
            subplot(2,2,3); imagesc(xz);  colorbar;  hold on; plot(fitResults.x,fitResults.z,'ro');   title('data xz');
            subplot(2,2,4); imagesc(xzS); colorbar;  hold on; plot(par(f+1:4:end),par(f+3:4:end),'ro');  title('fit xz');
            colormap(gray);
        else
            figure(pars.figShowFits); clf;
            [xy,xz,yz] = ProjectIm3D(data_3d); % use the original data
            subplot(3,1,1); imagesc(xy);  colorbar; hold on; plot(fitResults.x,fitResults.y,'ro'); title('data xy');
            subplot(3,1,2); imagesc(xz);  colorbar;  hold on; plot(fitResults.x,fitResults.z,'ro'); title('data xz');
            subplot(3,1,3); imagesc(yz);  colorbar;  hold on; plot(fitResults.y,fitResults.z,'ro'); title('data yz');
            colormap(gray);
        end
    end
    