
% Build Blast Libs
TSTORMlib = 'C:\Users\Alistair\Documents\Research\Projects\TSTORM\Lib\';
usedFlyPrimers = [TSTORMlib, 'UsedFlyPrimers.fasta'];
secondaries = [TSTORMlib, 'secondaries.fasta'];
usedTSTORMprimers = [TSTORMlib, 'UsedTSTORMprimers.fasta']; % includes secondaries
imr90transcriptome = [TSTORMlib, 'TotRNArep1totFPKM.fasta'];
rRNAandTRNA = [TSTORMlib, 'HumanRrnaAndTrna.fasta'];


% BuildBLASTlib(usedTSTORMprimers,'legacy',false);
% BuildBLASTlib(rRNAandTRNA,'legacy',false);
% BuildBLASTlib(imr90transcriptome,'legacy',false);



BLASTlibs = {usedTSTORMprimers,imr90transcriptome,rRNAandTRNA};
maxBLASThomology = [13,13,10]; % max homology per lib

saveName = 'C:\Users\Alistair\Documents\Research\Projects\TSTORM\Lib\Lib5primers.fasta';

BuildOrthogPrimers(saveName,'BLASTlibs',BLASTlibs,'maxBLASThomology',maxBLASThomology);