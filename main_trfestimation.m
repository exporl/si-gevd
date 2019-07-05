% MAIN_TRFESTIMATION Builds synthetic EEG and checks the quality of TRF estimation:
% Raw Vs CCA-filtered Vs SI-GEVD filtered
%
% This script builds synthetic EEG based on base temporal response function (TRF) templates.
% Half the data with Left attended TRF and half the data with right attended TRF.
% Noise is added to it, and TRFs are estimated from the noisy EEG based on forward modelling.
% These are compared with when TRFs are estimated filtering of the noisy EEG data with an SI-GEVD
% filter or a CCA filter. The results are compared by using at therelative MSE with the base template.
% Author: Neetha Das
% Date: 6/6/2019

clc;
clear all;

addpath(genpath(pwd))

% Set number of GEVD or CCA components to use
comps = 2;
%Sampling frequency
fs = 32;
%Trial length
tlen = 120;
no_samples = ceil(tlen*fs);
lambda = 0.2; % regularization parameter

% Set range of stimulus lag (in ms)
lag_start = 0;
lag_end = 400;

%Range of SNR
SNR_targets =  0:-5:-25;


% Load the base TRFs (Left and Right), and the EEG data+stimulus for a subject (15_subject in the AAD dataset) to be used as noise
load('data_TRFestim.mat');
eegnoise = flipud(subjectdata.eeg); % flip in time to ensure decorrelation with stimulus
env = subjectdata.stimenv;
assert(subjectdata.fs == fs);
TRFL = meanTRFL;
TRFR = meanTRFR;

% Find the channel amplitudes at the reference lag - at which the TRFs jointly have the largest norm across channels
for i = 1:14
    nrm(i) = norm([TRFL(i,:) TRFR(i,:)]);
end
[~,maxlag] = max(nrm);
amp_L_maxlag = abs(TRFL(maxlag,:));
amp_R_maxlag = abs(TRFR(maxlag,:));

% Use TRF of channel C5 as a base template
ch = 14;
base_TRF_L = TRFL(:,ch);
base_TRF_R = TRFR(:,ch);

% Build the full base TRFs by scaling the base template with the amplitude at the refence lag
base_TRF_L_all = repmat(base_TRF_L,1,63)*diag(amp_L_maxlag);
base_TRF_R_all = repmat(base_TRF_R,1,63)*diag(amp_R_maxlag);

%Use the audio envelope to predict stimulus following responses (EEG space)
idx_s = floor(lag_start/1e3*fs);
idx_f = ceil(lag_end/1e3*fs);
stim_lags = idx_s:idx_f;
lagged_stim = lag_data(env,stim_lags)';

EEGatt_L = base_TRF_L'*lagged_stim;
EEGatt_L_allch = repmat(EEGatt_L,63,1)';
EEGatt_L_allch = EEGatt_L_allch*diag(amp_L_maxlag);

EEGatt_R = base_TRF_R'*lagged_stim;
EEGatt_R_allch = repmat(EEGatt_R,63,1)';
EEGatt_R_allch = EEGatt_R_allch*diag(amp_R_maxlag);

% Split EEG into trials
signal_length = size(eegnoise,1);
no_splits = floor(signal_length/no_samples);
no_splits_all = no_splits*2;% X 2 because left attended as well as right attended TRFs are separately convolved with the same stimulus

mse_raw = zeros(size(SNR_targets,2),no_splits_all);
mse_gevd = zeros(size(SNR_targets,2),no_splits_all);
mse_cca = zeros(size(SNR_targets,2),no_splits_all);

figure(5);
ha = tight_subplot(2,3,[.05 .01],[.1 .05],[.05 .01]);

% For a range of SNRs simulate noisy EEG, and find TRFs on a per trial basis
for SNRnum = 1: size(SNR_targets,2)
    
    SNR_target =  SNR_targets(SNRnum);
    
    % scaling to achieve the target SNR
    SNR = 10^(SNR_target/20);
    rms_noise = rms(eegnoise);
    
    rms_EEG = rms(EEGatt_L_allch);
    SNR_base = mean(rms_EEG)/mean(rms_noise);
    C_L = SNR_base/SNR;
    
    rms_EEG = rms(EEGatt_R_allch);
    SNR_base = mean(rms_EEG)/mean(rms_noise);
    C_R = SNR_base/SNR;
    
    % Add noise to the synthesized EEG
    respL = EEGatt_L_allch + C_L*eegnoise;
    respR = EEGatt_R_allch + C_R*eegnoise;
    
    stim = env;
    
    count = 0;
    total_range = no_splits_all:-1:1;
    
    % Total data includes equal amount of 120s trials for L and R attended cases
    resp = [respL(1:no_splits*no_samples,:); respR(1:no_splits*no_samples,:)];
    stim = [stim(1:no_splits*no_samples,:); stim(1:no_splits*no_samples,:)];
    att_dir = [ones(no_splits*no_samples,1); 2*ones(no_splits*no_samples,1)]; % Direction of attention: 1 - Left, 2 - Right
    signal_length = size(resp,1);
    
    % Split EEG into trials and find GEVD and CCA filters based on a LOO basis, from the training set
    for j = total_range
        
        count = count + 1;
        range = signal_length-j*no_samples+1:signal_length-(j-1)*no_samples;
        assert(length(range)==ceil(tlen*fs));
        
        % Ensure that attention was to the same direction during the full trial
        assert(att_dir(range(1)) == att_dir(range(end)))
        
        %Raw EEG
        EegData = resp(range,:);
        
        %EEG for training the GEVD and CCA filters
        EegData_train = resp;
        EegData_train(range,:) = []; % leaving a trial out, to build the training set
        stim_train = stim;
        stim_train(range,:) = [];
        att_dir_train = att_dir;
        att_dir_train(range,:) = [];
        
        % Normalize data
        stim_train=normalize_data(stim_train);
        EegData_train=normalize_data(EegData_train);
        
        % Find the denoising filters
        [U,W,~] = find_sigevd_filter(EegData_train,stim_train,att_dir_train,fs,lag_start,lag_end,comps,lambda); %,no_of_loops);
        [Ucca,Wcca,stim_lags] = find_cca_filter(EegData_train,stim_train,fs,lag_start,lag_end,comps,lambda);
        
        % Test set
        resp_test = resp(range,:);
        stim_test = stim(range,:);
        
        %Apply the denoising filters
        EegData_gevd =  resp_test*W;
        EegData_cca =  resp_test*Wcca;
        
        %Audio Data
        lagged_stim = lag_data(stim_test,stim_lags);
        t = stim_lags(1)/fs:1/fs:stim_lags(end)/fs;
        
        %Estimate the TRFs from
        %Raw EEG
        trf_raw = trfestim(lagged_stim',EegData',lambda);
        %GEVD filtered EEG
        trf_gevd = trfestim(lagged_stim',EegData_gevd',lambda);
        %CCA filtered EEG
        trf_cca = trfestim(lagged_stim',EegData_cca',lambda);
        
        % Find the base TRF from which the test trial was simulated
        if att_dir(range(1)) == 1
            base_TRF_all = base_TRF_L_all;
        else
            base_TRF_all = base_TRF_R_all;
        end
        
        % Scale the estimated TRFs to the same range as the base TRF (using least squares)
        alpha_raw = trf_raw(:)\(base_TRF_all(:));
        trf_raw_scaled = trf_raw*alpha_raw;
        
        alpha_gevd = trf_gevd(:)\(base_TRF_all(:));
        trf_gevd_scaled = trf_gevd*alpha_gevd;
        
        alpha_cca = trf_cca(:)\(base_TRF_all(:));
        trf_cca_scaled = trf_cca*alpha_cca;
        
        % Find the relative mean squared error between the base TRF and the estimated TRF
        mse_raw(SNRnum,count) = (norm(base_TRF_all(:) - trf_raw_scaled(:))/norm(base_TRF_all(:)))^2;
        mse_gevd(SNRnum,count) = (norm(base_TRF_all(:) - trf_gevd_scaled(:))/norm(base_TRF_all(:)))^2;
        mse_cca(SNRnum,count) = (norm(base_TRF_all(:) - trf_cca_scaled(:))/norm(base_TRF_all(:)))^2;
        
        if j == 1 % for one trial per SNR, show the TRFs
            
            axes(ha(SNRnum));
            ha(SNRnum).XTickLabelMode = 'auto';
            t = t.*1000;
            line(t',base_TRF_all(:,ch),'LineWidth',1.5,'Color',[0 0 1],'Linestyle','-')
            hold on;
            line(t',trf_raw_scaled(:,ch),'LineWidth',1.5,'Color',[1 0 0],'Linestyle',':')
            line(t',trf_cca_scaled(:,ch),'LineWidth',1.5,'Color',[0.5 0 0.5],'Linestyle','--')
            line(t',trf_gevd_scaled(:,ch),'LineWidth',1.5,'Color',[0 0.75 0],'Linestyle','-.')
            xlim([0 410])
            if SNRnum > 3
                xlabel('Lags (ms)','FontSize',12)
            end
            if SNRnum == 3
                legend({'Base','Raw','CCA','SI-GEVD'},'FontSize',12)
            end
            title([num2str(SNR_target) 'dB'],'FontSize',12)
            
        end
    end
    
end

set(ha(1:3),'XTickLabel',''); set(ha,'YTickLabel','')

%% Move data to an MS Excel spreadsheet (statistical analysis and boxplots in the SI-GEVD paper were made in R)
if 1
    data1 = {};
    data2 = {};
    data3 = {};
    count = 1;
    
    for SNRnum = 1:length(SNR_targets)
        for j = 1: size(mse_raw,2)
            data1(count,:) = {num2str(j), num2str(tlen),num2str(SNR_targets(SNRnum)), 'SI-GEVD',num2str(mse_gevd(SNRnum,j))};
            data2(count,:) = {num2str(j), num2str(tlen),num2str(SNR_targets(SNRnum)), 'Raw',num2str(mse_raw(SNRnum,j))};
            data3(count,:) = {num2str(j), num2str(tlen),num2str(SNR_targets(SNRnum)), 'CCA',num2str(mse_cca(SNRnum,j))};
            count = count + 1;
        end
        
    end
    
    data_header = [{'Trialid','Triallength','SNR','Method','relMSE'};
        data1;
        data2;
        data3];
    xlswrite([pwd filesep 'relMSE_LR_TRFs_' num2str(tlen) 's.xlsx'],data_header);
    
end

