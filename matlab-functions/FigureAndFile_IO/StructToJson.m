function StructToJson(structIn,filename)
% saves structIn as filename, where structIn is a structure and filename is
% a the file path, which can be .txt or .json 
fid = fopen(filename,'w+');
jdata = jsonencode(structIn,'PrettyPrint',true)
fprintf(fid,jdata);
fclose(fid);