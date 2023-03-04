function structOut = JsonToStruct(jsonFile)

% load json and decode
fid = fopen(jsonFile,'r');
raw = char(fread(fid)');
fclose(fid);
structOut = jsondecode(raw);