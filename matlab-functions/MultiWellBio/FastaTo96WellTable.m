function FastaTo96WellTable(fastaIn,tableName)

if ischar(fastaIn)
    fastaIn = fastaread(fastaIn);
end

Name = {fastaIn.Header}';
Sequence = {fastaIn.Sequence}';
WellPosition = Index96Well();
WellPosition = WellPosition(1:length(Name));
T = table(WellPosition,Name,Sequence);
writetable(T,tableName); 