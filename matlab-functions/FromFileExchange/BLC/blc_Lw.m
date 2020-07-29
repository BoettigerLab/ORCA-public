function [LL, gradient, H] = blc_Lw(w,alpha,model)

% [LL, gradient, H] = blc_Lw(w,alpha,model)
%
% Initialises a model structure.
%
% Inputs:
%
%   w = d+1 by 1 vector where the first d components are regression weights
%   and the d+1 component is w0, the bais term.
%   alpha = scalar valued hyperparamter controlling variance of weights.
%   model = model structure.
%
% Output:
%
%   LL = negative log of posterior over weights.
%   gradient = d by 1 vector of partial deriviatives at w
%   H = Hessian matrix at w
%
%
% Copyright 2013 James Barrett
%
% Date: 21 November 2013
% Contact: james.j.barrett@kcl.ac.uk

APPROX = 20;
BETA = 10;     % Variance of w0 prior

% Extract w0
w0 = w(model.d+1);
w(model.d+1) = [];

% Define h as w*x + w0
hminus = (w'*model.Xminus')' + w0;
hplus = (w'*model.Xplus')' + w0;

% Generate indices of individuals who require the approximation
minus_approx = find(hminus>=APPROX);
minus_ok = find(hminus<APPROX);
plus_approx = find(hplus<=-APPROX);
plus_ok = find(hplus>-APPROX);
hsqminus = hminus(minus_approx).^2;
hsqplus = hplus(plus_approx).^2;


% Contribution from individuals in -1 class
Lminus = zeros(model.Nminus,1);
Lminus(minus_ok) = log(0.5*erfc(hminus(minus_ok)));
Lminus(minus_approx) = -(hminus(minus_approx).^2) - log(hminus(minus_approx)) - log(2*sqrt(pi));
Lminus = sum(Lminus);

% Contribution from individuals in +1 class
Lplus = zeros(model.N-model.Nminus,1);
Lplus(plus_ok) = log(0.5*erfc(-hplus(plus_ok)));    % Equivalent to log(1 - 0.5*erfc(hplus(plus_ok))) by more stable
Lplus(plus_approx) = -(hplus(plus_approx).^2) - log(-hplus(plus_approx)) - log(2*sqrt(pi)); %Note that -h is used for log
Lplus = sum(Lplus);

% Prior terms
Lprior = -(0.5/(alpha^2)) * (w'*w) -...
    (0.5/(BETA^2)) * (w0^2);

% Constant terms
%Lconstant = (0.5*model.N * log(2*pi)) + log(alpha)

LL = -(1/model.N) * (Lminus + Lplus + Lprior);

%% Compute Gradient
if (nargout > 1)
    
    gradient = zeros(1,model.d+1);
    
    % 1 by d terms
    minus_term = zeros(model.Nminus,1);
    minus_term(minus_ok) = -(2/sqrt(pi)) *( (exp(-hminus(minus_ok).^2)) ./ (erfc(hminus(minus_ok))) );
    minus_term(minus_approx) = -2*hminus(minus_approx).* (1 -(2*hsqminus).^(-1) + 2*(2*hsqminus).^(-2)).^(-1);
    
    plus_term = zeros(model.N-model.Nminus,1);
    plus_term(plus_ok) = (2/sqrt(pi)) * (exp(-hplus(plus_ok).^2)) ./ (erfc(-hplus(plus_ok)));
    plus_term(plus_approx) = -2*hplus(plus_approx).* (1 -(2*hsqplus).^(-1) + 2*(2*hsqplus).^(-2)).^(-1);
    
    gradient(1,1:model.d) = sum(bsxfun(@times,minus_term,model.Xminus)) +...  % Contributions from -1 class
        sum(bsxfun(@times,plus_term,model.Xplus));                            % Contributions from +1 class
    
    gradient(1:model.d) =  gradient(1:model.d) - (w'/(alpha^2));   % Prior terms
    
    % Gradient of w0
    gradient(model.d+1) = sum(minus_term) + sum(plus_term) - (w0/(BETA^2));
    
    gradient = -(1/model.N) * gradient';
    
end

%% Compute Hessian
if (nargout > 2)
    
    A = zeros(model.d+1,model.d+1);    
    Aminus_term = -((minus_term.^2) + 2*hminus.*minus_term);
    Aplus_term = -((plus_term.^2) + 2*hplus.*plus_term);
    
    % d by d w terms
    term1 = bsxfun(@times,Aminus_term,model.Xminus);
    term2 = bsxfun(@times,Aplus_term,model.Xplus);
    for mu = 1:model.d
        A(mu,1:model.d) = sum(bsxfun(@times,model.Xminus(:,mu),term1)) +...
            sum(bsxfun(@times,model.Xplus(:,mu),term2));
    end
    A(1:model.d,1:model.d) = A(1:model.d,1:model.d) -(1/alpha^2)*eye(model.d);     % Prior terms
    
    % w0 cross terms
    A(model.d+1,1:model.d) = sum(bsxfun(@times,Aminus_term,model.Xminus))+...
        sum(bsxfun(@times,Aplus_term,model.Xplus));
    A(1:model.d,model.d+1) = A(model.d+1,1:model.d);  % Since matrix is symmetric
    
    % w0 self term
    A(model.d+1,model.d+1) = sum(Aminus_term)+sum(Aplus_term) - (1/BETA^2);
    
    H = -(1/model.N) * A;
    
end