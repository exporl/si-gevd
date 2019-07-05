function [Pk,W,stim_lags] = find_cca_filter(resp,stim,fs,start,fin,no_of_comps,lambda)
% FIND_CCA_FILTER Estimate the CCA filter
%Input:
% resp: Neural response - time samples X channels
% stim: Stimulus - time samples X 1
% fs: Sampling frequency (Hz)
% start: start time lag of TRF in milliseconds
% fin: stop time lag of TRF in milliseconds
% no_of_comps: number of components from the CCA to be used
%Output:
% Pk: Spatial filters from EEG to CCA subspace - channels X components
% W: Denoising filter - channels X channels
% Author: Neetha Das
% Date: 6/6/2019

orig_resp = resp;

% lags for stimulus
idx_s = floor(start/1e3*fs);
idx_f = ceil(fin/1e3*fs);
stim_lags =idx_s:idx_f;
stim = lag_data(stim,stim_lags);

% Center the variables
[n,~] = size(stim);
stim = stim - repmat(mean(stim,1), n, 1);
resp = resp - repmat(mean(resp,1), n, 1);

S = stim'; % Lags X time
R = resp'; % Channels X time

XXT = S*S'; % stim covar matrix
XXT = XXT + eye(size(XXT)).*lambda*max(max(abs(XXT)));

invXXT = inv(XXT);
XYT = S*R';

YYT = R*R'; % EEG covariance matrix
invYYT = inv(YYT);
YXT = R*S';

% solution for the CCA filter for EEG
[B,r] = eig(invYYT*YXT*invXXT*XYT);
[~,desc_idx] = sort(real(diag(r)),'descend');
B = B(:,desc_idx);

% Regression to find the back projection filter
Pk = B(:,1:no_of_comps);
x = Pk'*resp'; % dim: components X time
% find Qk such that Qk'*x is as close as possible to resp'
Qk = (x*x')\(x*orig_resp); %inv(x*x')*x*resp;
W = Pk*Qk;

end