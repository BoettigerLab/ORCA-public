function indexNames = Index96Well
% Returns a 96x1 cell array containing the names A01-H12

indexNames = cell(96,1);
indexLetters = {'A','B','C','D','E','F','G','H'};
p = 0;
for l = 1:8
    for n=1:12
        p=p+1;
        indexNames{p} = [indexLetters{l},sprintf('%02d',n)];
    end
end