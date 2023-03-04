function fitTable = GaussFit2D(im2,dTable,varargin)

% parameters inherited from FitPsf2D
defaults = cell(0,3);
defaults(end+1,:) = {'cameraBackground','nonnegative',0};
defaults(end+1,:) = {'minPeakHeight','nonnegative',1};
defaults(end+1,:) = {'peakBlur','nonnegative',0};
defaults(end+1,:) = {'trimBorder','nonnegative',0};
defaults(end+1,:) = {'troubleshoot','boolean',false};
defaults(end+1,:) = {'cropBoxWidth','nonnegative',8};
defaults(end+1,:) = {'keepBrightest','integer',inf};
% fitter
defaults(end+1,:) = {'FiniteDifferenceType',{'central','forward'},'central'};
defaults(end+1,:) = {'MaxFunctionEvaluations','positive',5E3};
defaults(end+1,:) = {'OptimalityTolerance','positive',1E-8};
defaults(end+1,:) = {'method',{'fit','lsqn'},'fit'};
% units
defaults(end+1,:) = {'xyUnitConvert','positive',1};
pars = ParseVariableArguments(varargin,defaults,mfilename);


% handle vector input instead of table
if ~istable(dTable)
    xy = dTable;
    x = xy(:,1);
    y = xy(:,2);
    dTable = table(x,y);
end

% complete table 
if ~isfield(dTable,'peakHeight')
    idx = sub2ind(size(im2),dTable.x,dTable.y);
    dTable.peakHeight = im2(idx);
end

nPts = height(dTable);
x = zeros(nPts,1);
y = zeros(nPts,1);
wx = zeros(nPts,1);
wy = zeros(nPts,1); 
h = zeros(nPts,1); 
a = zeros(nPts,1); 
b = zeros(nPts,1); 
xL = zeros(nPts,1); 
xU = zeros(nPts,1); 
yL = zeros(nPts,1); 
yU = zeros(nPts,1); 
for i = 1:nPts
    cropWidth = pars.cropBoxWidth; 
    [v,w] = size(im2);
    xr = max(1,dTable.x(i)-cropWidth):min(dTable.x(i)+cropWidth,w);
    yr = max(1,dTable.y(i)-cropWidth):min(dTable.y(i)+cropWidth,v);
    data_2d = double(im2(yr,xr)); % 
   
    [rows,cols] = size(data_2d);
    [Y,X] = meshgrid(1:cols, 1:rows);

    a0 = double(dTable.peakHeight(i));
    b0 = min(data_2d(:));
    mu_x0 = dTable.x(i)-min(xr)+1; %  cropWidth+1;
    mu_y0 = dTable.y(i)-min(yr)+1; % cropWidth+1;
    sigma_x0 = 1;
    sigma_y0 = 1;
    peakBound = 2;
    minSigma = .1;
    maxSigma = 3;
 
    if strcmp(pars.method,'lsqn')
        %--------------------------------------------------------------------%  
        % least squares fitting to a 3D Gaussian
        minRes = @(p) exp(-((Y(:)-p(1))/(2*p(2))).^2 -((X(:)-p(3))/(2*p(4))).^2 )*p(5)+p(6) -data_2d(:);
        par0 = [mu_x0,sigma_x0,mu_y0,sigma_y0,a0,b0]; % initial fits
        lb = [mu_x0-peakBound,minSigma,mu_y0-peakBound,minSigma,-inf,-inf]; % lower bound 
        ub = [mu_x0+peakBound,maxSigma,mu_y0+peakBound,maxSigma,inf,inf]; % upper bound
        options = optimoptions('lsqnonlin',...
                               'FiniteDifferenceType',pars.FiniteDifferenceType,...
                               'OptimalityTolerance',pars.OptimalityTolerance,...
                               'MaxFunctionEvaluations',pars.MaxFunctionEvaluations,...
                               'Display','off');
        [par,~,residual,~,~,~,jacobian] = lsqnonlin(minRes, par0,lb,ub,options);   
        ci = nlparci(par,residual,'jacobian',jacobian,'alpha',1-.95); % 95% CI 

        x(i) = (min(xr)-1+par(1))*pars.xyUnitConvert;
        y(i) = (min(yr)-1+par(3))*pars.xyUnitConvert;
        wx(i) = par(2)*pars.xyUnitConvert;
        wy(i) = par(4)*pars.xyUnitConvert;
        h(i) = dTable.peakHeight(i);
        a(i) = par(5);
        b(i) = par(6); 
        xL(i) = (min(xr)-1+ci(1,1))*pars.xyUnitConvert;
        xU(i) = (min(xr)-1+ci(1,2))*pars.xyUnitConvert;
        yL(i) = (min(yr)-1+ci(3,1))*pars.xyUnitConvert;
        yU(i) = (min(yr)-1+ci(3,2))*pars.xyUnitConvert;
        % resRatio(i) = resnorm/a(i)^2;
    %---------------------------------------------------------------------%
    elseif strcmp(pars.method,'fit')
        ftype = fittype('exp(-((x-mu_x)/(2*sigma_x)).^2-((v-mu_y)/(2*sigma_y)).^2 )*a +b',...
                        'coeff', {'a','mu_x','sigma_x','mu_y','sigma_y','b'},'ind',{'x','v'});   % a = a/(2*pi*sigma_x*sigma_y)

        % dataTest = a0*exp(-((X-mu_x0-1.3)/(2*sigma_x0+3)).^2-((Y-mu_y0+.25)/(2*sigma_y0+.5)).^2 )+b0;
        fit_2d = fit([Y(:),X(:)],double(data_2d(:)),ftype,'StartPoint',[a0 mu_x0 sigma_x0 mu_y0 sigma_y0 b0],...
            'Lower',[-inf mu_x0-peakBound minSigma mu_y0-peakBound minSigma -inf],...
            'Upper',[ inf mu_x0+peakBound maxSigma mu_y0+peakBound maxSigma inf]);
        bounds = confint(fit_2d); 

        x(i) = min(xr)-1+fit_2d.mu_x;
        y(i) = min(yr)-1+fit_2d.mu_y;
        wx(i) = fit_2d.sigma_x;
        wy(i) = fit_2d.sigma_y;
        h(i) = dTable.peakHeight(i);
        a(i) = fit_2d.a;
        b(i) = fit_2d.b;
        xL(i) = bounds(1,2)+min(xr)-1;
        xU(i) = bounds(2,2)+min(xr)-1;
        yL(i) = bounds(1,4)+min(yr)-1;
        yU(i) = bounds(2,4)+min(yr)-1;  
    %-------------------------------------------------------------------%
    end
    
end


fitTable = table(x,y,wx,wy,h,a,b,xL,xU,yL,yU);



