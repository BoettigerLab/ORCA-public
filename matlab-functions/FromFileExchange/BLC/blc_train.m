function [w, w0, werror, w0error, logL] = blc_train(X,sigma)
% Bayesian linear classifier.
%
% [w, w0, werror, w0error, logL] = blc_train(X,sigma)
%
% Trains a bayesian liner classfier where
%
%               p(sigma=+1|x) = 0.5*erfc(w'*x + w0)
%
% and
%
%               p(sigma=-1|x) = 1 - p(sigma=+1|x)
%
% Inputs:
%
%   X = N by d matrix of covariates where N is the number of samples
%   (or observations) and d is the number of covairates (or input
%   variables).
%   sigma = N by 1 vector of class labels '+1' or '-1'
%
% Outputs:
%
%   w = d by 1 vector of regression weights
%   w0 = bias term
%   werror = estimate of standard deviation of w (based on curvature matrix
%       around minimum)
%   w0error = standard deviation of w0
%   logL = log likelihood
%
% See blc_demo.m for an example.
% See blc_documentation.pdf for further details.
%
%
% Copyright 2013 James Barrett
%
% Date: 21 November 2013
% Contact: james.j.barrett@kcl.ac.uk

% Set random number seed
s = RandStream('mcg16807','Seed',314159);
RandStream.setGlobalStream(s);

% Initialise model
model = blc_initialise(X,sigma);

% Set Options for Optimisation
options = optimset;
options.GradObj = 'on';
options.Hessian = 'on';
options.Display = 'off';
options.MaxIter = 10*model.N;
options.TolX = 1e-6;
options.TolFun = 1e-6;

% Optimise log likelihood with respect to w
initial_w = randn(model.d+1,1);
alpha = 1;
wstar = fminunc(@(w)blc_Lw(w,alpha,model), initial_w, options);
w0 = wstar(model.d+1);
w = wstar(1:model.d);

% Get error bars on w
if (nargout > 1)
[logL, ~, H] = blc_Lw(wstar,alpha,model);
L = chol(model.N*H);    % Use cholesky decomposition to invert A
invL = inv(L);
invA = invL*invL';
werror = sqrt(diag(invA));
w0error = werror(model.d+1);
werror(model.d+1)=[];
end