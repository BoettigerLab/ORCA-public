function [spots,pars] = ChrTracer3p2_ManualSelectSpot(im,varargin)

% defaults
defaults = cell(0,3);
defaults(end+1,:) = {'displayContrastLow','fraction',0.3};
defaults(end+1,:) = {'displayContrastHigh','fraction',.999};
defaults(end+1,:) = {'lociXY','freeType',[]};
defaults(end+1,:) = {'verbose','boolean',true};
pars = ParseVariableArguments(varargin,defaults,mfilename);

imOut = IncreaseContrast(im,'low',pars.displayContrastLow,'high',pars.displayContrastHigh);
overlayFig = figure(1); clf; 
imagesc(imOut); hold on;
colormap(gray); 

if ~isempty(pars.lociXY)
    plot(pars.lociXY(:,1),pars.lociXY(:,2),'yo');
    spotNums = cellstr(num2str((1:size(pars.lociXY,1))'));
    text(pars.lociXY(:,1)+1,pars.lociXY(:,2),spotNums,'color','w');
end

% prompt user to select locus of interest
dcmObj = datacursormode(overlayFig);
set(dcmObj,'DisplayStyle','datatip','SnapToDataVertex','off','Enable','on')
disp('Select a spot with the Data Cursor tool, then press Return.')
pause  % Wait while the user does this.hand
c_info = getCursorInfo(dcmObj);
if pars.verbose
    disp('Position captured'); 
end
spots = c_info.Position;