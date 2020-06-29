function [dTable,parameters] = FitPsf3D(Im1,varargin)
% Description:
% FitPsf3D takes the input matrix, finds all local maxima above a certain
% threshold, and then fits 3D Gaussians to these.  
% Someday this should be rewritten in C.  
% 
% -------------------------------------------------------------------------
% Required Inputs
% -------------------------------------------------------------------------
% Im1, a VxWxZ 3D-image
% 
% -------------------------------------------------------------------------
% Outputs
% -------------------------------------------------------------------------
% fits, a table containing the following fields
% fits.x - x position of peak
% fits.y - y position of peak
% fits.z - z position of peak
% fits.wx - sigma_x for the PSF
% fits.wy - sigma_y for the PSF
% fits.wz - sigma_z for the PSF
% fits.h - the pixel intensity of local region max pixel
% fits.a  - the fitted peak height of the unnormalized Gaussian
% fits.b - the fitted base height of the unnormalized Gaussian
% 
% -------------------------------------------------------------------------
% Default Options
% -------------------------------------------------------------------------
% % key parameters
% defaults(end+1,:) = {'minPeakHeight', 'positive', 1000};
% defaults(end+1,:) = {'cameraBackground', 'positive', 0};
% defaults(end+1,:) = {'peakBlur', 'positive', .5};
% defaults(end+1,:) = {'maxFitWidth', 'positive', 8};
% defaults(end+1,:) = {'initSigmaXY','positive',1.25};
% defaults(end+1,:) = {'initSigmaZ','positive',2.5};
% % less common parameters
% defaults(end+1,:) = {'minSep','positive',3}; % minimum separation in pixels between initial maxima
% defaults(end+1,:) = {'maxSigma','positive',2.5}; % max psf XY (PSF sigmaZ ~= 2*sigmaXY). Upperbound in fit
% defaults(end+1,:) = {'minSigma','positive',.1}; % min psf XY (PSF sigmaZ ~= 2*sigmaXY). Lower bound in fit
% defaults(end+1,:) = {'peakBound','positive',2}; % max number of pixels between brightest-pixel and gaussian-center. Set upper and lower bounds during fitting
% defaults(end+1,:) = {'FiniteDifferenceType',{'central','forward'},'central'};
% defaults(end+1,:) = {'MaxFunctionEvaluations','positive',5E3};
% defaults(end+1,:) = {'OptimalityTolerance','positive',1E-8};
% defaults(end+1,:) = {'verbose','boolean',true};
% defaults(end+1,:) = {'troubleshoot','boolean',true};
% 
% -------------------------------------------------------------------------
% Notes
% 
% 
% -------------------------------------------------------------------------
% Alistair Boettiger (boettiger@stanford.edu)
% February 14, 2017
% Copyright CC BY NC
% -------------------------------------------------------------------------
% 
% Updates
% 8/7/17  Split off 3D peak finding to a separate function, FindPeaks3D
% which is called by this function.  Left 3D Gaussian fitting intact.
% 
% 
%%

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
% key parameters
defaults(end+1,:) = {'minPeakHeight', 'positive', 1000};
defaults(end+1,:) = {'cameraBackground', 'nonnegative', 0};
defaults(end+1,:) = {'peakBlur', 'nonnegative', .5}; % gaussian smoothing before initial find max 
defaults(end+1,:) = {'maxFitWidth', 'positive', 8};
defaults(end+1,:) = {'initSigmaXY','positive',1.25};
defaults(end+1,:) = {'initSigmaZ','positive',2.5};
% less common parameters
defaults(end+1,:) = {'minSep','positive',3}; % minimum separation in pixels between initial maxima
defaults(end+1,:) = {'maxSigma','positive',2.5}; % max psf XY (PSF sigmaZ ~= 2*sigmaXY). Upperbound in fit
defaults(end+1,:) = {'minSigma','positive',.1}; % min psf XY (PSF sigmaZ ~= 2*sigmaXY). Lower bound in fit
defaults(end+1,:) = {'peakBound','positive',2}; % max number of pixels between brightest-pixel and gaussian-center. Set upper and lower bounds during fitting
defaults(end+1,:) = {'FiniteDifferenceType',{'central','forward'},'central'};
defaults(end+1,:) = {'MaxFunctionEvaluations','positive',5E3};
defaults(end+1,:) = {'OptimalityTolerance','positive',1E-8};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'veryVerbose','boolean',false};
defaults(end+1,:) = {'troubleshoot','boolean',true};
% filtering parameters
defaults(end+1,:) = {'keepBrightest','integer',inf};
defaults(end+1,:) = {'relativeHeight','fraction',0};
defaults(end+1,:) = {'minHBratio','nonnegative',1.2}; % peak value over background value
defaults(end+1,:) = {'minAHratio','nonnegative',.25}; % fitted height over background vs peak value
defaults(end+1,:) = {'maxUncert','nonnegative',2}; % pixels
defaults(end+1,:) = {'resRatio','nonnegative',inf}; % pixels
% manual seed point
defaults(end+1,:) = {'seedPoint','array',[]}; 
% units
defaults(end+1,:) = {'xyUnitConvert','positive',1};
defaults(end+1,:) = {'zUnitConvert','positive',1};

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

warning('off','MATLAB:singularMatrix');

% -------------------------------------------------------------------------
% Main Function
% -------------------------------------------------------------------------

% Find peaks
if isempty(parameters.seedPoint) 
   dTable = FindPeaks3D(Im1,'parameters',parameters,'troubleshoot',false);
else
    x = parameters.seedPoint(1);
    y = parameters.seedPoint(2);
    z = parameters.seedPoint(3);
    h = Im1(x,y,z); 
    dTable = table(x,y,z,h);
end
w = floor(parameters.maxFitWidth/2); % short hand for crop;

if parameters.troubleshoot
    % plot local maxima to pixel precision
    subplot(3,2,1); imagesc(max(Im1,[],3)); hold on; plot(dTable.x,dTable.y,'o','color',[1 .5 .5]);  % plots x vs y
    % squashed the columns (y-dim); put the stack (z-dim) on the y-axis (new columns);  plots x vs z;
    subplot(3,2,2); imagesc(max(permute(Im1,[3,2,1]),[],3)); hold on; plot(dTable.x,dTable.z,'o','color',[1 .5 .5]); 
end

% Fit 3D Gaussian
nMax = length(dTable.x);
wx = zeros(nMax,1);
wy = zeros(nMax,1);
wz = zeros(nMax,1);
a = zeros(nMax,1);
b = zeros(nMax,1); 
xL = zeros(nMax,1); 
xU = zeros(nMax,1); 
yL = zeros(nMax,1); 
yU = zeros(nMax,1); 
zL = zeros(nMax,1); 
zU = zeros(nMax,1); 
resRatio = zeros(nMax,1); 

for i = 1:nMax
    % crop region around local maxima
    [rows,cols,stcks] = size(Im1);
    xr = max(1,dTable.x(i)-w):min(dTable.x(i)+w,cols);
    yr = max(1,dTable.y(i)-w):min(dTable.y(i)+w,rows);
    zr = max(1,dTable.z(i)-w):min(dTable.z(i)+w,stcks);
    data_3d = double(Im1(yr,xr,zr));  % matrix indexing, col-y, row-x, stck-z
    
    if parameters.troubleshoot
        % plot region being used for spot-fitting
        [Y,X] = meshgrid(yr,xr);  % matrix indexing, col-y, row-x, stck-z
        subplot(3,2,3); imagesc(max(Im1,[],3)); hold on; plot(X(:),Y(:),'r.');  % plot indexing, x vs y 
        subplot(3,2,4); imagesc(max(data_3d,[],3)); 
        title(i);
    end
    
    % parameters for fit
    % [cols,rows,stcks] = size(data_3d);
    [rows,cols,stcks] = size(data_3d);
    [Y,X,Z] = meshgrid(1:cols, 1:rows, 1:stcks);
    a0 = double(dTable.h(i));
    b0 = 300; % quantile(data_3d(:),.1);
    mu_x0 = dTable.x(i)-min(xr)+1;
    mu_y0 = dTable.y(i)-min(yr)+1; 
    mu_z0 = dTable.z(i)-min(zr)+1;
    sigma_x0 = parameters.initSigmaXY;
    sigma_y0 = parameters.initSigmaXY;
    sigma_z0 = parameters.initSigmaZ;
    maxSigma = parameters.maxSigma; 
    minSigma = parameters.minSigma;
    peakBound = parameters.peakBound;

    % least squares fitting to a 3D Gaussian
    minRes = @(p) exp(-((Y(:)-p(1))/(2*p(2))).^2 -((X(:)-p(3))/(2*p(4))).^2 -((Z(:)-p(5))/(2*p(6))).^2 )*p(7)+p(8) -data_3d(:);
    par0 = [mu_x0,sigma_x0,mu_y0,sigma_y0,mu_z0,sigma_z0,a0,b0]; % initial fits
    lb = [mu_x0-peakBound,minSigma,mu_y0-peakBound,minSigma,mu_z0-peakBound,minSigma,0,0]; % lower bound 
    ub = [mu_x0+peakBound,maxSigma,mu_y0+peakBound,maxSigma,mu_z0+peakBound,2*maxSigma,2^16,2^16]; % upper bound
    options = optimoptions('lsqnonlin',...
                           'FiniteDifferenceType',parameters.FiniteDifferenceType,...
                           'OptimalityTolerance',parameters.OptimalityTolerance,...
                           'MaxFunctionEvaluations',parameters.MaxFunctionEvaluations,...
                           'Display','off');
    [par,resnorm,residual,~,~,~,jacobian] = lsqnonlin(minRes, par0,lb,ub,options);
    
    ci = nlparci(par,residual,'jacobian',jacobian,'alpha',1-.95); % 95% CI 
    
    % resnorm is just: resnorm = sum( (fit(:) - data(:)).^2 )
    % can also get residuals = fit(:) - data(:).  
    % what we really want are confidence interverals on the mu_x,y,z.
    
    % save data
    dTable.x(i) = (min(xr)-1+par(1))*parameters.xyUnitConvert;
    dTable.y(i) = (min(yr)-1+par(3))*parameters.xyUnitConvert;
    dTable.z(i) = (min(zr)-1+par(5))*parameters.zUnitConvert;
    wx(i) = par(2)*parameters.xyUnitConvert;
    wy(i) = par(4)*parameters.xyUnitConvert;
    wz(i) = par(6)*parameters.zUnitConvert;
    a(i) = par(7);
    b(i) = par(8); 
    xL(i) = (min(xr)-1+ci(1,1))*parameters.xyUnitConvert;
    xU(i) = (min(xr)-1+ci(1,2))*parameters.xyUnitConvert;
    yL(i) = (min(yr)-1+ci(3,1))*parameters.xyUnitConvert;
    yU(i) = (min(yr)-1+ci(3,2))*parameters.xyUnitConvert;
    zL(i) = (min(zr)-1+ci(5,1))*parameters.zUnitConvert;
    zU(i) = (min(zr)-1+ci(5,2))*parameters.zUnitConvert;
    resRatio(i) = resnorm/a(i)^2;
    
    % plot results for validation
    if parameters.troubleshoot
        subplot(3,2,5); imagesc(max(data_3d,[],3));  
        hold on; plot(mu_x0,mu_y0,'ro'); plot(par(1),par(3),'r+');
        subplot(3,2,6); imagesc(max(permute(data_3d,[3,2,1]),[],3));  
        hold on; plot(mu_x0,mu_z0,'ro'); plot(par(1),par(5),'r+');
    end
end


dTable = [dTable,table(wx,wy,wz,a,b,xL,xU,yL,yU,zL,zU)];

% filter 'spotlike' parameters
keepSpots = (dTable.xU - dTable.xL < parameters.maxUncert*parameters.xyUnitConvert) ...
            & (dTable.yU - dTable.yL < parameters.maxUncert*parameters.xyUnitConvert) ...
            & (dTable.zU - dTable.zL < 2*parameters.maxUncert*parameters.zUnitConvert);
keepSpots = keepSpots & dTable.h ./ dTable.b >= parameters.minHBratio;
keepSpots = keepSpots & dTable.a ./ dTable.h >= parameters.minAHratio;
keepSpots = keepSpots & resRatio < parameters.resRatio;


if parameters.troubleshoot
    dTable %#ok<NOPRT>
    resRatio %#ok<NOPRT>
    keepSpots  %#ok<NOPRT>   
end
dTable = dTable(keepSpots,:);


if parameters.troubleshoot  && nMax > 0
    subplot(3,2,1);  hold on; plot(dTable.x/parameters.xyUnitConvert,dTable.y/parameters.xyUnitConvert,'r+');
    plot(min(xr)-1+mu_x0,min(yr)-1+mu_y0,'m+'); 
    subplot(3,2,2);  hold on; plot(dTable.x/parameters.xyUnitConvert,dTable.z/parameters.zUnitConvert,'r+');
    plot(min(xr)-1+mu_x0,min(zr)-1+mu_z0,'m+'); 
end

