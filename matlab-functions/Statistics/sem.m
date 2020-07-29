function s = sem(v)
% return the standard error of the mean

s = nanstd(v)/sqrt(length(v));