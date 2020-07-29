function y = FlatProtectedDecay(t, A, a, k)
%--------------------------------------------------------------------------
% y = FlatProtectedDecay(t, A, a, k)
% This function returns amount of mRNA at time t subjected to the protected
% decay model with amplitude A, transcription time a, and decay rate k.
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% January 27, 2015
%
% Version 1.0
%--------------------------------------------------------------------------
% Creative Commons License CC BY
%--------------------------------------------------------------------------
ind = t<a;
y(ind) = 1;
ind = t>=a;
y(ind) = exp(-k*t(ind))*exp(k*a);

y = A*y';
