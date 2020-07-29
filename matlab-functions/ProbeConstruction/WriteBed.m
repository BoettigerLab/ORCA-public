function datafile = WriteBed(bedTable,saveFolder,fname,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'genome','string',''};
defaults(end+1,:) = {'description','string',''};
defaults(end+1,:) = {'browserPosition','string',''};
defaults(end+1,:) = {'removeTable','boolean',true};
defaults(end+1,:) = {'verbose','boolean',true};
pars = ParseVariableArguments(varargin,defaults,mfilename);

if isempty(pars.description)
    description = [pars.genome,' ',fname];
else
    description = pars.description;
end

header = {};
if ~isempty(pars.browserPosition)
   header{end+1} = ['browser position: ',pars.browserPosition]; 
end
header{end+1} = ['track name="',fname,'" description="',description,'" visibility=1 itemRgb="On"'];
bedfile =[saveFolder,fname,'_bed.txt']; 
headerfile = [saveFolder,fname,'_header.txt'];
datafile = [saveFolder,fname,'_data.txt'];
if exist(bedfile,'file')~=0
    delete(bedfile);
end

fid = fopen(headerfile,'w+');
fprintf(fid,'%s\r\n',header{:});
fclose(fid); 
writetable(bedTable,datafile,'delimiter','\t','WriteVariableNames',false);

[~,~] = system(['copy /b ',headerfile,'+',datafile,' ',bedfile]);
if pars.verbose
    disp(['wrote ',bedfile]);
end
if pars.removeTable
    delete(datafile);
end
delete(headerfile);