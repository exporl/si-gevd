function [lagged_data] = lag_data(data,lags)
% LAG_DATA Adds lagged channels (with zero padding) to the data matrix
%Inputs:
% data: data matrix (time x channels)
% lags: time-lags to be used
%Outputs:
% lagged_data: lag data matrix (time x (channels x lags))
% Author: Neetha Das
% Date: 6/6/2019

dim = size(data,2);
lagged_data = zeros(size(data,1),dim*length(lags));
idx = 1;

for i = 1:length(lags)
    
    % Shift data
    lag_mat = circshift(data,lags(i));
    
    % Zero-pad
    if lags(i) < 0
        lag_mat(end-abs(lags(i))+1:end,:) = 0;
    else
        lag_mat(1:lags(i),:) = 0;
    end
    
    % Update output matrix for this lag
    lagged_data(:,idx:idx+dim-1) = lag_mat(1:size(data,1),:);
    idx = idx + dim;
end

end