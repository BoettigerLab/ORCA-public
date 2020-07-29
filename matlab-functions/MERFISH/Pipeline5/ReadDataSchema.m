function headerStruct = ReadDataSchema(experimentLayout)
% just reads headerData
% 
% Alistair Boettiger
% 
%%
% experimentLayout = '\\morgan\MorganData2\MERFISHdata\150709_L16Test3\10uM\dataSchema_150721_L16_10uM.xlsx';

headerData = readtable(experimentLayout,'ReadVariableNames',false,'ReadRowNames',true,'Range','A1:B5');
fieldNames = regexprep(headerData.Properties.RowNames,' ',''); % make names legal  
headerStruct = cell2struct(headerData.Var1,fieldNames,1);

