function [p,h_ratio] = SurvivalComp(a,b,varargin)
% Compare the means of two assumed exponential distributions 
% 
% coerce the Cox proportional hazards model into the same format as
% ranksum, kstest2, ttest2, etc. 
% 
% https://en.wikipedia.org/wiki/Proportional_hazards_model#References
% 
% a = cd_all_dTag{7}';
% b = cd_all_dTag{8}';
a = Column(a);
b = Column(b);
X = [0*ones(length(a),1); 1*ones(length(b),1)];
T = [a; b];
if ~isempty(varargin)
mdl = fitcox(X,T,'Censoring',varargin{1});
else
mdl = fitcox(X,T);
end
p= mdl.Coefficients.pValue;
h_ratio = exp(mdl.Coefficients.Beta);
