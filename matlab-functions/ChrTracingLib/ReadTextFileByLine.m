function text = ReadTextFileByLine(txtFile)    
% returns a cell array containing the data on each line of txtFile as a
% string

fid = fopen(txtFile);
count = 1;
text = {};
while ~feof(fid)
    text{count} = fgetl(fid);
    count = count + 1;
end
fclose(fid);