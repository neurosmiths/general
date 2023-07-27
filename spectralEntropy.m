function [Hs,nBins,theoreticalMax] = spectralEntropy(X)
% SPECTRALENTROPY normalizes and calculates the entropy of a frequency spectrum
%
%   [dE] = spectralEntropy(spectrum) calculates the spectral entropy
%   of the spectrum input, after normalizing the spectrum to PSD.  
%	
%	[dE,nBins,theoreticalMax] = spectralEntropy(signal) will also output the 
%	number of bins used for spectral estimation and the theoretical maximum
%	value: 
%		
%		theoreticalMax = log2(nBins);
% 

% getting number of bins in spectrum
nBins = length(X);

% normalizing spectrum
spectralDist = X./sum(X);

% calculating entropy
Hs = -nansum(spectralDist.*log2(spectralDist));

% calculating the theoretical maximum. 
theoreticalMax = log2(nBins);
