function [rho,p] = LogCorr(x,y)

skip = ~(isnan(x) | isnan(y) | x==0 | y==0 | isinf(x) | isinf(y));
[rho,p] = corr(log10(x(skip)),log10(y(skip)));