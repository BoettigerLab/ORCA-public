function datPropTable = DataChnsFromTable(eTable)

numHybes = height(eTable);
channels = cellfun(@(x) strsplit(regexprep(x,'[^A-Za-z0-9,]',''),','),eTable.channels,'UniformOutput',false);
fidChannel = strcmp(channels{1},num2str(eTable.fiducialChannel(1)));
datChannels = channels{1}(~fidChannel);
numDatChns = length(channels{1})-1;
numDats = numHybes*numDatChns;
readout = zeros(numDats,1);
chn = cell(numDats,1);
dataType = cell(numDats,1);
hybe = zeros(numDats,1);
chnNum = zeros(numDats,1);
k=0;
try
for h=1:numHybes
    if iscell(eTable.Readouts) % new format, a string of multiple readouts, comma separated
        readNames = strsplit(eTable.Readouts{h},','); 
    else % old format, a number
        readNames = string(eTable.Readouts(h)); 
    end
    for r=1:numDatChns
        k=k+1;
        readout(k) = str2num(readNames{r}); %#ok<ST2NM>
        chn{k} = datChannels{r};
        dataType{k} = eTable.DataType{h};
        hybe(k) = h;
        chnNum(k) = r;
    end
end
datPropTable = table(readout,chn,dataType,hybe,chnNum);
catch er
    disp(er.getReport);
    disp('debug here');
end