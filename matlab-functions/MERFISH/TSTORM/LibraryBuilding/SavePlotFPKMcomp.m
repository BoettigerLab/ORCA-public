function SavePlotFPKMcomp(saveName,data1,varargin)

if nargin > 2
    data2 = varargin{1}; 
end


figure(1); clf; colordef white; set(gcf,'color','w');
loglog(data1.FPKM,data1.FPKM,'k.');  hold on;
loglog(data2.FPKM,5*data2.FPKM,'r.');
for i=1:length(data1.CommonName)
text(data1.FPKM(i)+.15*data1.FPKM(i),(data1.FPKM(i)),data1.CommonName{i},'color','k');
end
for i=1:length(data2.CommonName)
text(data2.FPKM(i)+.15*data2.FPKM(i),(5*data2.FPKM(i)),data2.CommonName{i},'color','r');
end
xlabel('FPKM');
export_fig(saveName); 