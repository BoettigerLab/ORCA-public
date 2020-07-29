function PlotExpectHybe(libExpectPerHybe,hybeLocs)

numHybes = length(hybeLocs); 
plot([libExpectPerHybe*mean(hybeLocs)*numHybes,hybeLocs],'.-')
PresentationPlot('MarkerWidth',25);
xlabel('hybe');
ylabel('# localizations'); 
legend('Expected','Observed');
xlim([1,numHybes]);