function order = ShortestPath3D(points,varargin)


defaults = cell(0,3);
defaults(end+1,:) = {'verbose','boolean',true};  %
defaults(end+1,:) = {'christofidesPath','boolean','U:\Alistair\MonetBackup\Research\Software\Lemon_mathematica\'};  %
pars = ParseVariableArguments(defaults,varargin,mfilename);

[N,dims] = size(points);
if dims == 2
    points = [points,zeros(N,1)];
end
fprintf('Create/Load points\n');
christofidesPath = pars.christofidesPath; 
christofides = [christofidesPath,'christofides'];

%% Save positions to file
filename=[christofidesPath,'test.dat'];
fileID = fopen(filename,'w');
fprintf(fileID,'%s\n',filename);
fprintf(fileID,'%i\n',N);

for i=1:N
    fprintf(fileID,'%d\n%d\n%d\n',points(i,1),points(i,2),points(i,3));
end
fclose(fileID);

%% Call christofides.exe with data
tic;
[~,log]=system([christofides,' < ',filename]);
% system([christofides,' < ',filename])  % I think this is redundant
toc;
C=strsplit(log,'\n');
C=strsplit(C{end-1},{',',':'});
order=str2double({C{2:end-1}})+1;
