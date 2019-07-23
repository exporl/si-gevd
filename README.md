
# SI-GEVD

## License

See the LICENSE file for license rights and limitations. By downloading and/or installing this software and associated files on your computing system you agree to use the software under the terms and condition as specified in the License agreement.

## The SI-GEVD implementation

### About

This MATLAB code implements an algorithm based on the SI-GEVD filtering for denoising multi-channel EEG as presented in the Validation Experiment `Short-term TRF estimation’ in [1]. SI-GEVD filtering improves the signal-to-noise ratio (SNR) of the stimulus following responses in the EEG data. Stimulus-following EEG data at different SNRs (0 dB to -25 dB) are first simulated. From 120 s trials, temporal response functions (TRFs) are estimated for 3 cases: from raw data, from SI-GEVD-filtered data, and from canonical correlation analysis (CCA)-filtered data. The relative mean square errors (relMSEs) with respect to the base TRF templates are computed as a measure of the quality of TRF estimation in each case.
Developed and tested in MATLAB R2015a. 

### Documentation

All functions are documented properly in their respective m-files. The main file is main_trfestimation.m, which calls the other functions in order to implement the algorithm and display and save the results. The boxplots of relMSEs as shown in [1] can be generated from plot_results.R. Additional information about how the EEG data was synthesized, how SI-GEVD filter was designed, and the similarities and differences with respect to CCA etc. can be found in [1].

### References

[1] Das, N., Vanthornhout, J., Francart, T. and Bertrand, A., 2019. Stimulus-aware spatial filtering for single-trial neural response and temporal response function estimation in high-density EEG with applications in auditory research. bioRxiv, p.541318.
