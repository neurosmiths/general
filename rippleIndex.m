function [I] = rippleIndex(X,f)
% RIPPLEINDEX normalizes and calculates the fast ripple index of a spectrum.
%
%   [I] = rippleIndex(X) calculates the fast ripple index
%   of the spectrum in vector X, sample at the frequencies in f. 
%
% 	This is based on the calculation described in Ibarz et al. (2010) JNeurosci.
%

if ~isequal(length(f),length(X))
	error('X and f must have the same number of elements')
end

% normalizing spectrum
spectralDist = X./sum(X);

% calculating RI
I = trapz(X(f>90))./trapz(X);




