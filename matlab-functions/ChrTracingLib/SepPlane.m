function [crossA,crossB,estimates] = SepPlane(setA,setB,varargin)
% [crossA,crossB,estimates] = SepPlane(setA,setB,varargin)
% compute the best separation plane between the points setA and setB and
% return a logical vector of the points of A which lie on the side of B and 
% the points from B which lie on the side dominated by A (crossA, crossB). 
% optionally also return the plane itself (estimates, the a,b,c,d values)
% from the plane equation 0 = a*X+b*Y+c*Z+d
% 
% % Example
% nPts = 100;
% setA = rand(nPts,3);  
% setB = 0.1 + rand(nPts,3) + repmat([.3,.1,.2],nPts,1);

defaults = cell(0,3);
defaults(end+1,:) = {'showPlots','boolean',false};
pars = ParseVariableArguments(varargin,defaults,mfilename);

% deal with NaNs/missing data
skipA = isnan(setA(:,1));
setA(skipA,:) = [];
skipB = isnan(setB(:,1));
setB(skipB,:) = [];

% ---------- compute separation plane -----------
estimates = SeparationPlane(setA, setB);

% -------- ID points on either side ------------
% Point-to-Plane signed distance approach
% distance from plane
% if point = (x1,y1,z1) 
% and plane: 0 = a*X+b*Y+c*Z+d
% distance point to plane = (a*x1+b*y1+c*z1+d) / sqrt(a^2+b^2+c^2)
a = estimates(1); b= estimates(2); c=estimates(3); d=estimates(4);
setA_dist = (a*setA(:,1)+b*setA(:,2)+c*setA(:,3)+d) / sqrt(a^2+b^2+c^2);
setB_dist = (a*setB(:,1)+b*setB(:,2)+c*setB(:,3)+d) / sqrt(a^2+b^2+c^2);
crossA = setA_dist > 0;
crossB = setB_dist < 0;

% ------------- show plots for troubleshooting -------------
if pars.showPlots
    figure(1); clf; 
    plot3(setA(:,1),setA(:,2),setA(:,3),'r.'); hold on;
    plot3(setB(:,1),setB(:,2),setB(:,3),'b.');
    xlabel('x'); ylabel('y'); zlabel('z');

    xmin = min([setA(:,1); setB(:,1)]);
    xmax = max([setA(:,1); setB(:,1)]);
    ymin = min([setA(:,2); setB(:,2)]);
    ymax = max([setA(:,2); setB(:,2)]);
    [X,Y] = meshgrid(linspace(xmin,xmax,100),linspace(ymin,ymax,100));
    % ax + by + cz + ... + d = 0
    Z = -(estimates(1)*X+estimates(2)*Y+ estimates(4))/estimates(3);
    hold on;
    surf(X,Y,Z);
    
    plot3(setA(crossA,1),setA(crossA,2),setA(crossA,3),'mo');
    plot3(setB(crossB,1),setB(crossB,2),setB(crossB,3),'go');
    title( sum(~crossA)/length(crossA) );
end

% restore nan's
crossA_out = nan(length(skipA),1);
crossA_out(~skipA) = crossA;
crossA = crossA_out;
crossB_out = nan(length(skipB),1);
crossB_out(~skipB) = crossB;
crossB = crossB_out;