function[trf] = trfestim(lagged_stim, resp, lambda)
% TRFESTIM Estimate temporal response function
%Input:
% lagged_stim: Stimulus matrix - lags X time
% resp: Neural response - channels X time
% lambda: regularisation parameter
%Output:
% trf: Temporal response function - lags X channels
% Author: Neetha Das
% Date: 6/6/2019

XXT = lagged_stim*lagged_stim';
XYT = lagged_stim*resp';

% Set up regularisation
XXT = XXT + eye(size(XXT)).*lambda*max(max(abs(XXT)));

% Find the TRF
trf = XXT\XYT;

end