function [Hd,theoreticalMax] = differentialEntropy(X,n)
% DIFFERENTIALENTROPY calculates the differential entropy of a signal
%
%   [dE] = differentialEntropy(signal) calculates the differential entropy
%   of the continuous signal in vector X.  
%
%	[dE,theoreticalMax] = differentialEntropy(signal,nBins) also outputs the 
%	theoretical maximum value for differential entropy, given the number
%	of bins used, nBins, The default number of bins is 100. 

% default number of bins
if ~exist('n','var')
    n = 100;
end

% the random distirbution for comparison
k = 1; % fortifies the uniform distribution with extra data.
mX = histcounts(rand(1,length(X)*k),n)./(length(X)*k);

% calculating differential entropy. 
[p,edges] = histcounts(X,n);
p = p./length(X);
Hd = -nansum(p.*log2(p));

% calculating theoretical maximum value
theoreticalMax = log2(n);

