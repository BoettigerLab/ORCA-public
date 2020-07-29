function y = ProtectedDecay(t, A, a, k)
%--------------------------------------------------------------------------
% y = ProtectedDecay(t, A, a, k)
% This function returns amount of mRNA at time t subjected to the protected
% decay model with amplitude A, transcription time a, and decay rate k.
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% February 27, 2013
%
% Version 1.0
%--------------------------------------------------------------------------
% Creative Commons License CC BY
%--------------------------------------------------------------------------
ind = t<a;
y(ind) = a*(1-(t(ind)/a).^2) + 1/k;
ind = t>=a;
y(ind) = 1/k*exp(-k*t(ind))*exp(k*a);

y = A*y';
