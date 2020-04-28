function sym = Symmetrize(sym)
% symmetrize a matrix
% short functions written originally for aparna's nueral net paper
% 
% 

sym = triu(sym)+tril(sym)';
sym = (sym + sym')/2;
