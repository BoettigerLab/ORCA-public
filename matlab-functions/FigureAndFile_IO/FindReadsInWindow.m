function regTable = FindReadsInWindow(locustxt,dataTable,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'spacer','nonnegative',0};

parameters = ParseVariableArguments(varargin,defaults,mfilename);

% locustxt = 'chr3R:12481406-12810708'; % BXC
% dataTable = table(K79me3_data{1},K79me3_data{2},K79me3_data{3},K79me3_data{4},...
%             'VariableNames',{'chr','locusStart','locusEnd','score'});

[chrName,locusStart,locusEnd] = ParseLocusName(locustxt);
locusStart = locusStart - parameters.spacer;
locusEnd = locusEnd + parameters.spacer;
chrId = StringFind(dataTable.chr,chrName);
chrTable = dataTable(chrId,:);

[err,startId] = min(abs(chrTable.locusStart - locusStart));
[err,endId] = min(abs(chrTable.locusEnd - locusEnd));

regTable = chrTable(startId:endId,:);

%  figure(2); clf; plot( mean([regTable.locusStart,regTable.locusEnd],2),...
%                         regTable.score,'k.');