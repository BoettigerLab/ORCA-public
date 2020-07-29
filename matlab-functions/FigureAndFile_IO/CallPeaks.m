function [peakLocs,peakHeights] = CallPeaks(x,y,varargin)
% [peakLocs,peakHeights] = CallPeaks(x,y)
% [peakLocs,peakHeights] = CallPeaks(x,y,'minHeight',30,'minDistance',30E3)
% 
% OPTIONS
% defaults(end+1,:) = {'minHeight','float',1};
% defaults(end+1,:) = {'showplots','boolean',true};
% defaults(end+1,:) = {'minDistance','nonnegative',0};

defaults = cell(0,3);
defaults(end+1,:) = {'minHeight','float',0};
defaults(end+1,:) = {'showplots','boolean',false};
defaults(end+1,:) = {'minDistance','nonnegative',0};
pars = ParseVariableArguments(varargin,defaults,mfilename);

keepData = y > pars.minHeight;
y = y(keepData);
x = x(keepData);
yy = smooth(x,y);  % probably should pass some smooth options 

[peakHeights,peakLocs] = findpeaks(yy,x,'MinPeakDistance',pars.minDistance); 

if pars.showplots
    plot(x,yy);
    hold on; plot(peakLocs,peakHeights,'v');
end