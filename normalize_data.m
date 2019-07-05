function [data] = normalize_data(data)
% NORMALIZE_DATA Set data to zero mean and unit standard deviation per channel
%Input:
% data: Stimulus/neural response - time X channels
%Output:
% data: Stimulus/neural response - time X channels
% Author: Neetha Das
% Date: 6/6/2019

data = data-repmat(mean(data),size(data,1),1);
stdev = std(data);
stdev(stdev<10e-5) = 10e-5;
data = data./repmat(stdev,size(data,1),1);

end