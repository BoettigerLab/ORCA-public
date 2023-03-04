function [chrNames,vals] = LoadWig(wigFile,varargin)

% wigFile =  'U:\GenomeData\Fly\Kc167\GSM1318352_Rad21.wig'
fid = fopen(wigFile);
fmt =['%f %f %*[^\n]'];
n=0;
headers = cell(100,1);
vals = cell(100,1);
while n<100
    n=n+1;
    header = textscan(fid,'%s',1,'delimiter','\t');
    disp(header{1})
    header = textscan(fid,'%s',1,'delimiter','\t');
    disp(header{1})
    headers{n} = header{1};
    val = textscan(fid,fmt,'CollectOutput',true,'delimiter','\t','HeaderLines',2); 
    vals{n} = val{1};
end
fclose(fid);
drop = cellfun(@isempty,headers);
headers(drop) = [];
vals(drop) = [];
datNames = cat(1,headers{:})
e = cellfun(@(x) strfind(x,'span'),datNames,'UniformOutput',false);
chrNames = cellfun(@(x,y) x(20:y(1)-2),datNames,e,'UniformOutput',false);

if ~isempty(varargin)
   locusTxt = varargin{1};
   [trk,xs] = ParseWig(vals,locusTxt,chrNames); 
   chrNames = trk; % overwrite output
   vals = xs; % overwrite output
end