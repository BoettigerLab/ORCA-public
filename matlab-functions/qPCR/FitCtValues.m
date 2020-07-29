function [cValues, cCI,fitResults] = FitCtValues(data)
% data is a matrix numCycles by numSamples

% data = standards;
[numCycles, numSamples] = size(data);

ft = fittype( 'a + b/(1+(2*d)^(-(x-c)))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( ft );
opts.Lower = [0 0 0 0];
opts.Upper = [Inf Inf Inf Inf];
xData = (1:numCycles)';
fits = {};
clear gofs
aValues = [];
bValues = [];
cValues = [];
dValues = [];
aCI = [];
bCI = [];
cCI = [];
dCI = [];

for i=1:numSamples
    yData = data(:,i);
    opts.StartPoint = [yData(1) yData(end) 15 1];
    [fitresult, gof] = fit( xData, yData, ft, opts );
    fits{i} = fitresult;
    aValues(i) = fitresult.a;
    bValues(i) = fitresult.b;
    cValues(i) = fitresult.c;
    dValues(i) = fitresult.d;
    ci = confint(fitresult);
    aCI(i,:) = ci(:,1)';
    bCI(i,:) = ci(:,2)';
    cCI(i,:) = ci(:,3)';
    dCI(i,:) = ci(:,4)';
    gofs(i) =gof;
end


% figure(5); clf;
% subplot(2,2,1); imagesc(reshape(cValues,12,8)); title('cValues'); colorbar;
% subplot(2,2,2); imagesc(reshape(aValues,12,8)); title('aValues'); colorbar; 
% subplot(2,2,3); imagesc(reshape(bValues,12,8)); title('bValues'); colorbar;
% subplot(2,2,4); imagesc(reshape(dValues,12,8)); title('dValues'); colorbar;
% 
% 
% figure(5); clf;
% subplot(2,2,1); imagesc(reshape(log10(cCI(:,2)),12,8)); title('cValues'); colorbar;
% subplot(2,2,2); imagesc(reshape(aValues,12,8)); title('aValues'); colorbar; 
% subplot(2,2,3); imagesc(reshape(bValues,12,8)); title('bValues'); colorbar;
% subplot(2,2,4); imagesc(reshape(dValues,12,8)); title('dValues'); colorbar;

reject = Column( (bValues - aValues)./aValues < 100 );
reject = reject | (cCI(:,2)-cCI(:,1))>3;
reject = reject | Column( cValues < 2 | cValues > numCycles-1 );

fitResults.aValues = aValues;
fitResults.bValues = bValues;
fitResults.cValues = cValues;
fitResults.dValues = dValues;

fitResults.aCI = aCI;
fitResults.bCI = bCI;
fitResults.cCI = cCI;
fitResults.dCI = dCI;


cValues(reject) = numCycles;

