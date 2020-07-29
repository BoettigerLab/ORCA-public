function model = blc_initialise(X,sigma)

% model = blc_initialise(X,sigma)
%
% Initialises a model structure.
%
%   X = N by d matrix of covariates where N is the number of samples
%   (or observations) and d is the number of covairates (or input
%   variables).
%   sigma = N by 1 vector of class labels '+1' and '-1'
%
% Output:
%
%   model = model structure
%
%
% Copyright 2013 James Barrett
%
% Date: 21 November 2013
% Contact: james.j.barrett@kcl.ac.uk

if (~isvector(sigma))
        error('blc:not_vector','class labels must be a column vector')
end

if (isrow(sigma))
    error('blc:not_column_vector','class labels must be a column vector')
end

checksigma = sigma;
checksigma(sigma==-1) = 1;
if (~all(checksigma==1))
    error('blc:class_labels','class labels must contain +1 or -1 only')
end

model.N = length(sigma);    %no. of individuals

if (abs(sum(sigma)) == model.N)
    error('blc:only_one_class','class labels cannot contain only one class')
end

if(size(X,1) ~= model.N)
    error('blc:dim_mismatch','input data must have the same number of rows as class labels')
end 

model.d = size(X,2);    %no. of covariates
model.X = X;
model.sigma = sigma;

%% Divide data into two classes
model.Nminus = sum(sigma==-1);
data = sortrows([X sigma],model.d+1);
model.Xminus = data(1:model.Nminus, 1:model.d);
model.Xplus = data(model.Nminus+1:model.N,1:model.d);
