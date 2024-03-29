% Script: cyclo_fdm_test.m
%
% Script showing cyclic spectrum for multiple FSK signal test file (fdm_test.mat). 
% User enters baud rate for analysis. Generates single cyclic spectrum
% (Cyclic Correlation Density (SCD)) from single baud rate input. 
%
% To Do: add correct calculations for adjusting noise to get desired SNR.
% 
% Author: drohm
%------------------------------------------------------------------------
%------------------------------------------------------------------------
clear all; close all;

%-Read in fdm test signal mat file
load fdm_test.mat;
fsk_c = x - mean(x);
td = [0:length(fsk_c)-1]*1/fs;

%-Cyclic frequency to test
alpha = input(['Input baud rate for analysis: ']);

%-Split out I and Q components of fsk signal
fsk_I = real(fsk_c);
fsk_Q = imag(fsk_c);

%-Generate AWGN noise
sigma_w = 100;                           % noise variance from input
Iwn = sigma_w*randn(1,length(fsk_I));
Qwn = sigma_w*randn(1,length(fsk_Q));
wn = Iwn + sqrt(-1)*Qwn;                % complex white gaussian noise

%-Add AWGN to complex baseband FSK signal
sigout = fsk_c + wn;

%-Generate Plots (optional)
plotswitch=1;
if plotswitch == 1
    figure
    
    subplot(2,1,1)
    plot(td,real(sigout))
    xlim([0 td(end)]);grid
    xlabel('Time')
    ylabel('Amplitude')
    title('FDM Waveform with AWGN')
   
    %-Take FFT of FSK signal
    y = sigout;
    nfft = 2.^(ceil(log(length(y))/log(2)));
    ffty = fft(y,nfft);                     
    ffty_mag = fftshift(abs(ffty));
    freq_limit = fs/2;
    freq= -freq_limit:fs/nfft:freq_limit-fs/nfft;
    
    %-Plot frequency domain
    subplot(2,1,2)
    % protect against taking log of zero or negative numbers
    dB_psd = 10*log10( (1/max(ffty_mag)) * ffty_mag + 1.e-10);
    plot(freq,dB_psd);    
    xlim([-freq_limit,freq_limit]);ylim([-50 10]);grid
    xlabel('FREQUENCY(Hz)');ylabel('DB');
    title('FSK Signal Spectrum')  
end
     
%-Compute frequency smoothed cyclic spectrum for selected baud rate
M = 8;                         % frequency smoothing factor
plotswitch = 1;
x = sigout(1:end);              % input signals (if x & y then compute cross-scd)
%[scd,freqx] = cyclodetect_fs(fs,alpha,M,plotswitch,x,y);

%-Compute time-smoothed cyclic spectrum for selected baud rate
[scd,freqx] = cyclodetect_ts(fs,alpha,M,plotswitch,x,y);


