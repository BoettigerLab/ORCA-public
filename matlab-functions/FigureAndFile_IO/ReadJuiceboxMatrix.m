function [plotMat,coords,max_coords] = ReadJuiceboxMatrix(matTxt,varargin)
% plotMat = ReadJuiceboxMatrix(matTxt,varargin)
%
%
% Example:
% matTxt = 'E:\Alistair\Analysis\2017-10-12_K562vsIMR90\IMR90_chr21reg_5kb_rawMatrix.txt';
% [plotMat,coords] = ReadJuiceboxMatrix(matTxt,...
%                    'locus','chr21:29372319-31372318');
% 
%

defaults = cell(0,3);
defaults(end+1,:) = {'showplot','boolean',false};
defaults(end+1,:) = {'locus','string',''};
defaults(end+1,:) = {'mapRes','nonnegative',5E3};
defaults(end+1,:) = {'minReads','nonnegative',0};
defaults(end+1,:) = {'displayRes','nonnegative',30E3};
pars = ParseVariableArguments(varargin,defaults,mfilename);



rawMat = readtable(matTxt);
if isempty(rawMat)
    error('table is empty! check filename');
end


stp =mode(diff(rawMat.Var1));
c1 = min( [rawMat.Var1; rawMat.Var2]);
c2 = max( [rawMat.Var1; rawMat.Var2]);
bins = (c2-c1)/stp+1;
x = (rawMat.Var1 - c1)/stp + 1;
y = (rawMat.Var2 - c1)/stp + 1;
max_coords = [c1,c2];

plotMat = zeros(bins,bins);
ind1 = sub2ind([bins,bins],x,y);
plotMat(ind1) = rawMat.Var3; % add UpperTriangle
ind2 = sub2ind([bins,bins],y,x);
plotMat(ind2) = rawMat.Var3; % add LowerTriangle
cm = quantile(rawMat.Var3,.95);

if ~isempty(pars.locus)
   [~,locusStart,locusEnd] = ParseLocusName(pars.locus);
   coords = [locusStart,locusEnd];
   c1new = ceil((locusStart - c1)/pars.mapRes);
   c2new = floor((locusEnd - c1)/pars.mapRes) ;
   if c1new < 1 || c2new < 1
      disp(['map region: ',num2str(c1),'-',num2str(c2),' request region: ',num2str(locusStart),'-',num2str(locusEnd)]);
      error(['locus ',pars.locus,' is not contained inside bounds of ', matTxt]); 
   end
   plotMat = plotMat(c1new:c2new,c1new:c2new);
else
    coords = max_coords;
end

plotMat = double(plotMat);

if pars.displayRes ~= pars.mapRes
    plotMat(plotMat==0) = nan;
    plotMat = InterpMapNans(plotMat);
    plotMat = plotMat./(max(plotMat(:)));
    sc = double(pars.mapRes/pars.displayRes);
    plotMat = double(imresize(uint16(2^15*plotMat),sc,'bicubic')); %    
end
if pars.showplot
    imagesc(log10(double(plotMat))); colorbar;
end