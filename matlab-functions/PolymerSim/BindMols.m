function Bound = BindMols(ka,kd,Bound,TotMols,varargin)

openChromatin = [];

% openChromatin,ka,FlagType,H3K27me3,kd)

N = length(Bound);
if isempty(openChromatin)
    openChromatin = ones(N,1); 
end
% ---------------Bind Proteins from Cytoplasm
freeMols = TotMols - sum(Bound); 
newBinders =  rand(N,1) < ka.*freeMols.*openChromatin*(~Bound); 
if sum(newBinders) > freeMols
  idx = datasample(find(newBinders),1,'Replace',false);
  newBinders = zeros(N,1); 
  newBinders(idx) = true;
end
Bound = Bound | newBinders;


% ------------ Release bound molecules at rate kd
unBound = rand(N,1) < kd.*(Bound); 
Bound = Bound & ~unBound; 
