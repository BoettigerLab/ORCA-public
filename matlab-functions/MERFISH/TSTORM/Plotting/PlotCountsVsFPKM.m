function pearsonCorr = PlotCountsVsFPKM(libExpect,EcCount,varargin)


x = libExpect; y = EcCount;
if ~isempty(varargin);
    pearsonCorr = PlotCorr(x,y,varargin{1});
else
    pearsonCorr = PlotCorr(x,y);
end
PresentationPlot('MarkerWidth',20);
xlabel('FPKM'); ylabel('counts');