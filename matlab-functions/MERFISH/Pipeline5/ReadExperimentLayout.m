function [bitOrder,bitNames,experimentInfo] = ReadExperimentLayout(experimentLayout)
% [idxBits,bitNames] = ReadExperimentLayout(experimentLayout)
% 
% Alistair Boettiger
% July 07, 2015

% experimentLayout = '\\Morgan\MorganData2\MERFISHdata\150703_L15L16_Test2\L15_10uM\ExperimentLayout.txt';

fid = fopen(experimentLayout);
experimentInfo = textscan(fid,'%s %s',1,'delimiter','\t');
fscanf(fid,'\n');
bitNameInfo = textscan(fid,'%s %s',1,'delimiter','\t');
fclose(fid);
bitNames = regexprep(bitNameInfo{2}{1},{'"','''',','},{'','',''});
bitNames = strsplit(bitNames,' ')';

bitTable = readtable(experimentLayout,'HeaderLines',3,'delimiter','\t');
idxBits = bitTable{:,2:end}; 
idxBits(isnan(idxBits)) = inf;
[numFrames,numHybes] = size(idxBits);
bitLabels = cell(numFrames,numHybes);
for f=1:numFrames
    for h=1:numHybes
        try
            if ~ischar(bitTable{f,1})
                bitName = num2str(bitTable{f,1});
            else
                bitName = bitTable{f,1};
            end
            bitLabels{f,h} = [bitNames{idxBits(f,h)},'-',bitName];
        catch
        end
    end
end

idxBits = idxBits(:);
bitNames = bitLabels(:);

[~,bitOrder] = sort(idxBits);
bitNames = bitNames(bitOrder);
