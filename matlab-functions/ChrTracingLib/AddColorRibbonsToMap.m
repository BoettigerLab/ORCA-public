function AddColorRibbonsToMap(im1,varargin)


defaults = cell(0,3);
defaults(end+1,:) = {'colormap', 'colormap', 'hsvCut'};
defaults(end+1,:) = {'MarkerSize', 'positive', 14};
pars = ParseVariableArguments(varargin,defaults,mfilename);

nBins = size(im1,1);
cmap = GetColorMap(pars.colormap,nBins); hold on;
for i=1:nBins
    plot([-.5,i],[i,-.5],'.','color',cmap(i,:),'MarkerSize',pars.MarkerSize);
end
xlim([-1,nBins]); ylim([-1,nBins]);

