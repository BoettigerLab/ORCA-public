function fovNums = GetFOVnums(daxNames1)

nFOV = length(daxNames1);
fovNums = zeros(1,nFOV);
for f=1:nFOV
    s = strfind(daxNames1{f},'_');
    s = s(end); % find last use;
    fnum = str2double(daxNames1{f}(s+1:end-4));
    fnum = fnum + 1 ; % convert from python 0 index
    fovNums(f) =  fnum; 
end
