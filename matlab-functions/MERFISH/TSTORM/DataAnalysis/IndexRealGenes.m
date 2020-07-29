function idxRealGenes = IndexRealGenes(libGenes)

blanks = StringFind(libGenes,'blank');
nontargets= StringFind(libGenes,'notarget');
idxRealGenes = 1:length(libGenes);
idxRealGenes([blanks; nontargets]) = [];
