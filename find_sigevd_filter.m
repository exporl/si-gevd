function [Pk,W,stim_lags] = find_sigevd_filter(resp,stim,att_dir,fs,start,fin,no_of_comps,lambda)
% FIND_SIGEVD_FILTER Estimate the SIGEVD filter
%Input:
% resp: Neural response - time samples X channels
% stim: Stimulus - time samples X 1
% att_dir: Direction of attention - time samples X 1 : values: 1 for Left, 2 for Right
% fs: Sampling frequency
% start: start time lag of TRF in milliseconds
% fin: stop time lag of TRF in milliseconds
% no_of_comps: number of components from the GEVD to be used
%Output:
% Pk: Generalized eigenvectors - channels X components
% W: Denoising filter - channels X channels
% Author: Neetha Das
% Date: 6/6/2019

% Split training data into L and R attended
idxL = find(att_dir == 1);
idxR = find(att_dir == 2);
respL = resp(idxL,:);
respR = resp(idxR,:);

stimL = stim(idxL,:);
stimR = stim(idxR,:);

% Apply lags for the stimulus
idx_s = floor(start/1e3*fs);
idx_f = ceil(fin/1e3*fs);
stim_lags = idx_s:idx_f;
stimL = lag_data(stimL,stim_lags);
stimR = lag_data(stimR,stim_lags);

% Center the variables
[n,~] = size(stimL);
stimL = stimL - repmat(mean(stimL,1), n, 1);
respL = respL - repmat(mean(respL,1), n, 1);

[n,~] = size(stimR);
stimR = stimR - repmat(mean(stimR,1), n, 1);
respR = respR - repmat(mean(respR,1), n, 1);

[n,~] = size(resp);
resp = resp - repmat(mean(resp,1), n, 1);

stimL = stimL';
stimR = stimR';
respL = respL';
respR = respR';
resp = resp';

% Estimate TRFs
trfa_L = trfestim(stimL, respL, lambda);
trfa_R = trfestim(stimR, respR, lambda);

% Stimulus following responses
Xa_L = trfa_L'*stimL;
Xa_R = trfa_R'*stimR;

% Signal covariance matrix
Ryy=(resp*resp')./size(resp,2);

% Stimulus following response convariance matrices - separately computed for left and right attended conditions
Ra_L = (Xa_L*Xa_L')./size(Xa_L,2);
Ra_R = (Xa_R*Xa_R')./size(Xa_R,2);
Ra = (Ra_L+Ra_R)./2;

% Perform GEVD
[U,lamb]=eig(Ra,Ryy);
[~,desc_idx] = sort(real(diag(lamb)),'descend');
U = U(:,desc_idx);

% Matrix for projection to GEVD subspace
Pk = U(:,1:no_of_comps);
Q = inv(U)';
% back projection filter
Qk = Q(:,1:no_of_comps);
% the full denoising filter
W = Pk*Qk';

end