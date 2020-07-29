function alignedIm = PlotAlignedDax(mRNAbin,listType,tform,varargin)
% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'showPlots', 'boolean', true};
defaults(end+1,:) = {'fighandle', 'handle', []};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 3
    error('matlabSTORM:invalidArguments', 'requires mRNAbin,listType,tform');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);
        

%% Main Function
numHybes = length(mRNAbin); 
snapshots = parameters.fighandle;  
showPlots = parameters.showPlots;

if isempty(snapshots) && showPlots
    snapshots = figure();
end
    
 
rawIm  = zeros(256,256,numHybes);
alignedIm = zeros(256,256,numHybes);
for h=1:numHybes
    daxfile = regexprep(mRNAbin{h},['_',listType],'.dax');
    dax = max(ReadDax(daxfile,'endFrame',1),[],3); % max project the first frame
    [H,W] = size(dax); 
    tformInv = fliptform(tform{h}); 
    alignedDax = imtransform(dax,tformInv,...
                    'XYScale',1,'XData',[1 W],'YData',[1 H]);
    rawIm(:,:,h) = dax;
    alignedIm(:,:,h) = alignedDax;
    if parameters.showPlots
        figure(snapshots); subplot(4,round(numHybes/4),h); 
        imagesc(uint16(alignedDax)); colormap gray;  % caxis([0,2^13.5]);
    end
end