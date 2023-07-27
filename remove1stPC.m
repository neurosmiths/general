function [denoisedData] = remove1stPC(data)
%REMOVE1STPC reconstructs data without the first PC.
%   common average re-referencing and de-meaning data 
%       (e.g. channels or trials by samples) with principal component analysis 
% 

% author: Elliot H Smith - https://github.com/elliothsmith/seizureCodes

% PCA
[w,pc,~] = pca(data);

% reconstruction without 1st component
denoisedData = (pc(:,2:end)*w(:,2:end)');

end

