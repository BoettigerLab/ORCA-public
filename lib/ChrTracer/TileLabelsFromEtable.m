function [tileLabels_fid,tileLabels_dat] = TileLabelsFromEtable(eTable)

% parse tile-labels from eTable
try
if ~isempty(eTable)
    if iscell(eTable.Bit(1))
        channels = cellfun(@(x) strsplit(regexprep(x,'[^A-Za-z0-9,]',''),','),eTable.channels,'UniformOutput',false);
        bitTypes = eTable.DataType;
        tileLabels_fid = cellfun(@(x,y) cat(2,x,y),eTable.Bit,bitTypes,'UniformOutput',false);
        temp = cellfun(@(x) strsplit(x,','),eTable.Bit,'UniformOutput',false);
        temp = cellfun(@(x,y,z) [strcat(x(1),'-',z(1),y), strcat(x(2),'-',z(2),y)],temp,bitTypes,channels,'UniformOutput',false);
        temp = cat(1,temp{:})';
        tileLabels_dat = temp(:);
        [numHybes,numDataChns] = size(tileLabels_dat);
        tileLabels_dat= reshape(tileLabels_dat,numHybes*numDataChns,1);
    else
        bitNums = cellstr(num2str(eTable.Bit));
        bitTypes = eTable.DataType;
        tileLabels = cellfun(@(x,y) cat(2,x,y),bitNums,bitTypes,'UniformOutput',false);
        tileLabels_fid = tileLabels;
        tileLabels_dat = tileLabels;
    end
else
    tileLabels_fid = [];
    tileLabels_dat = [];
end
catch
    tileLabels_fid = [];
    tileLabels_dat = [];
end