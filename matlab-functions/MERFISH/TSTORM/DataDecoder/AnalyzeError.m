function [ErrorPerHyb,FalseNegative,FalsePositive] = AnalyzeError(PmCount,codebook,uniqueMsg,hammingDis,freqUniqueMsg)
Ngenes = length(codebook);

Nbits = size(uniqueMsg,2);

for i=1:Ngenes
    missedbits(i) = sum(freqUniqueMsg((hammingDis(:,i)==1 & sum(uniqueMsg,2)==3)));
    gainedbits(i) = sum(freqUniqueMsg((hammingDis(:,i)==1 & sum(uniqueMsg,2)==5)));
end

codes = zeros(Ngenes,length(str2num(codebook(1).Header))); %#ok<*ST2NM>
for n = 1:Ngenes
    codes(n,:) = str2num(codebook(n).Header);
end 

% wash dependence of errors
washerr = zeros(Ngenes,Nbits);
for i=1:Ngenes
    onebiterrors =  uniqueMsg(hammingDis(:,i)==1,: );
    f = freqUniqueMsg(hammingDis(:,i)==1 );
    trueword = repmat(codes(i,:),size(onebiterrors,1),1);
    washerr(i,:) = f'*abs(trueword - onebiterrors);
end
for i=1:Nbits
    clusperhyb(i) = sum(uniqueMsg(:,i).*freqUniqueMsg);
end   
ErrorPerHyb = sum(washerr)./(clusperhyb+sum(washerr));
% figure(2); clf; plot(sum(washerr));
syms p
fn = 1-solve(p^4/(p^3*(1-p)*4) == (sum(PmCount)+sum(gainedbits))/sum(missedbits));
fp = 1-solve(p^4/(p^3*(1-p)*4) == (sum(PmCount)+sum(missedbits))/sum(gainedbits));
FalseNegative = double(fn);
FalsePositive = double(fp);

plot(ErrorPerHyb);
title(['false negative rate = ' num2str(FalseNegative) ', false positive rate = ' num2str(FalsePositive)])
xlabel('hybridization round');
PresentationPlot();

end