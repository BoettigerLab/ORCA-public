function [tileLabels_fid,tileLabels_dat] = TileLabelsFromEtable(eTable,varargin)
%%
% eTable = readtable('M:\2019-06-04_Rad21_TM3TwiGFP_18-20hrs_(HoxRNA)_(BXC3kbEven)\RNA_Expt\rnaTable_3color.xlsx')   
defaults = cell(0,3);
defaults(end+1,:) = {'style',{'readout','names'},'readout'};
pars = ParseVariableArguments(varargin,defaults,mfilename);

% parse tile-labels from eTable
try
if ~isempty(eTable)
    eTableS = table2struct(eTable);
    if strcmp(pars.style,'readout')    
        if isfield(eTableS,'Bit')
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
        elseif isfield(eTableS,'Readouts')
             if iscell(eTable.Readouts(1))         
                bitTypes = eTable.DataType;
                bitNum = cellstr(num2str(eTable.HybNum,'%02d'));
                tileLabels_fid =cellfun(@(x,y) strcat(x,'-',y),bitTypes,bitNum,'UniformOutput',false); %; cellfun(@(x,y) cat(2,x,y),eTable.Readouts,bitTypes,'UniformOutput',false);
                datPropTable = DataChnsFromTable(eTable);     
                tileLabels_dat =strcat(num2str(datPropTable.readout),datPropTable.dataType,'-',datPropTable.chn);     
            else
                tileLabels_fid = {};
                tileLabels_dat = {};
             end       
        end
    elseif strcmp(pars.style,'names')
        bitTypes = eTable.DataType;
        bitNum = cellstr(num2str(eTable.HybNum,'%02d'));
        tileLabels_fid =cellfun(@(x,y) strcat(x,'-',y),bitTypes,bitNum,'UniformOutput',false); %; cellfun(@(x,y) cat(2,x,y),eTable.Readouts,bitTypes,'UniformOutput',false);
              
        nameColumn = contains(eTable.Properties.VariableNames,'Name');
        nameColumn = find(nameColumn,1,'last');   
        chnNames = cellfun(@(x) strsplit(x,', '), eTable{:,nameColumn},'UniformOutput',false);
        chnNames = cellfun(@(x) cat(1,x(:)), chnNames,'UniformOutput',false);
        chnNames = cat(1,chnNames{:});
        tileLabels_dat = chnNames;
    end
else
    tileLabels_fid = {};
    tileLabels_dat = {};
end
catch
    tileLabels_fid = {};
    tileLabels_dat = {};
end