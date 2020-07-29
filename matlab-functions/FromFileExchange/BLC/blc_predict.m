function p = blc_predict(X,w,w0)

% p = blc_predict(X,w,w0)
%
% Calculates probability that a novel input X belongs to class +1
%
% Inputs:
%
%   X = M by d matrix of samples where each row is a sample
%   w = regression weights (typically from a trained model).
%   w0 = bias term.
%
% Output:
%
%   p = M by 1 vector of probabilities that each sample in X belong
%       to class +1. Computed according to
%
%                   sigma(i) = 1- 0.5*erfc(w'*X(:,i) + w0)
%
%
% Copyright 2013 James Barrett
%
% Date: 21 November 2013
% Contact: james.j.barrett@kcl.ac.uk

p = 1 - 0.5*erfc((w'*X')' + w0);