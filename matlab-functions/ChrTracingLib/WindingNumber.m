function [w,a] = WindingNumber(xy,varargin)
%% inputs
% xy - an Nx2 or Nx3 polymer array to compute the winding number of
%    - if 3D, x by y by z.
%% optional inputs
%  xy0, a 1x2 or 1x3 -- coordinates around which to compute the winding. 
% 
% Notes / source code
% https://comp.soft-sys.matlab.narkive.com/UDPLAsTQ/logicals-winding-number-inpolgon
% see reply 2 by Roger Stafford

if isempty(varargin)
    xy0 = [0,0];
else
    xy0 = varargin{1};
end

%% main funciton

% convert the inputs 
% xy = fillmissing(xy,'linear','EndValues','nearest');
xv = xy(:,1); yv=xy(:,2);
x0=xy0(1); y0=xy0(2);

xv = [xv;xv(1)]; yv = [yv;yv(1)]; % circularize (?)
x1 = xv(1:end-1); y1 = yv(1:end-1); % first coordinate
x2 = xv(2:end); y2 = yv(2:end);  % next coordinate 
a = atan2((x1-x0).*(y2-y0)-(y1-y0).*(x2-x0), ...
 (x1-x0).*(x2-x0)+(y1-y0).*(y2-y0))/(2*pi);
w = nansum(a);

%%  example
% theta = 0:.2:6*pi;
% x1 = 500*cos(theta)';
% y1 = 400*sin(theta)';
% 
% % wind the other way
% theta = 0:.2:4*pi;
% x2 = 300*cos(-theta)';
% y2 = 500*sin(-theta)';
% 
% x = [x1;x2]; y=[y1;y2];
% 
% 
% xy = [x,y];
% xy = fillmissing(xy,'linear','EndValues','nearest');
% [w,a] = WindingNumber([x,y])
% 
% figure(9); clf; bar(a); pause(.1);
% 
% figure(8); clf; 
% nB = length(xy); cmap = hsv(nB);
% plot(0,0,'k+'); hold on;xlim([-1000,1000]); ylim([-1000,1000]); 
% for i=1:nB
%     plot(xy(i,1),xy(i,2),'.','color', cmap(i,:)); 
%     if i<nB
%         plot(xy(i:i+1,1),xy(i:i+1,2),'.-','color', cmap(i,:));
%     end
%     title(nansum(a(1:i)))
%     pause(.3);
% end

