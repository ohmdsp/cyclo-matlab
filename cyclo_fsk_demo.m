% Script: cyclo_fsk.m
%
% Script showing cyclic spectrum for FSK signaling. Uses Matlab's 
% Communications Toolbox or FSK modulation.
% 
% Author: drohm
%------------------------------------------------------------------------
%------------------------------------------------------------------------
clear all; close all;

%-Generate a simulated FSK signal with AWGN
M = 2;                  % Modulation order
freqsep = 16;           % Frequency separation (Hz) (match to symbol rate)
nsamp = 8;             % Number of samples per symbol
fs = 128;               % Sample rate (Hz)
fcarr = 0;
rb = freqsep;
x = randi([0 M-1],1000,1);    % 200 symbols of max(td) is symbol rate
fsk_c = fskmod(x,M,freqsep,nsamp,fs,'cont'); 
td = [0:length(fsk_c)-1]*1/(fs);

%-Split out I and Q components of fsk signal
fsk_I = real(fsk_c);
fsk_Q = imag(fsk_c);

%-Generate AWGN noise
sigma_w = 0.02;                       % noise variance from input
Iwn = sigma_w*randn(1,length(fsk_I));
Qwn = sigma_w*randn(1,length(fsk_Q));
wn = Iwn + sqrt(-1)*Qwn;                % complex white gaussian noise

%-Add AWGN to complex baseband FSK signal
fsk_s = fsk_c + wn';

%-Generate carrier wave
twopi_fc_t = td*2*pi*fcarr;     % argument for carrier sinusoid 
a = 1;                          % amplitude
phi = 0;                        % initial phase
cs_t = a*exp(1i*twopi_fc_t + phi); % carrier signal
cs_t = cs_t';                   % transpose to column vector

%-Mix I/Q fsk data stream with carrier
sigout = fsk_s.*cs_t;           % mix symbols with carrier
sigout = sigout';            % reduce gain to keep output at +/- one

%-Generate Plots (optional)
plotswitch=1;
if plotswitch == 1
    figure
    
    subplot(2,1,1)
    plot(td,real(sigout))
    ylim([-1.5 1.5]); grid
    xlim([0 td(end)]);
    xlabel('Time')
    ylabel('Amplitude')
    title('FSK Waveform with AWGN')
   
    %-Take FFT of FSK signal
    y = sigout;
    nfft = 2.^(ceil(log(length(y))/log(2)));
    ffty = fft(y,nfft);             % take fft
    ffty_mag = abs(fftshift(ffty));
    freq_limit = fs/2;
    freq= -freq_limit:fs/nfft:freq_limit-fs/nfft;
    
    %-Plot frequency domain
    subplot(2,1,2)
    % this next line protects against taking log of zero or negative numbers
    dB_psd = 10*log10( (1/max(ffty_mag)) * ffty_mag + 1.e-10);
    plot(freq,dB_psd);
    xlim([-freq_limit,freq_limit]);ylim([-50 10]);grid
    xlabel('FREQUENCY(Hz)');ylabel('DB');
    title('FSK Signal Spectrum')  
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Cyclic Correlation    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-Compute non-conjugate cyclic autocorrelation of BPSK signal (2-D)
y = sigout;               % input signal (complex-valued)
max_lag = 32;             % number of +/- time lags in correlation
plotswitch = 1;           % generate plots 1->plots on, 0->plots off
r_alpha2d = ralpha_2d(y,fcarr,rb,max_lag,fs,plotswitch);

% %-Compute cyclic autocorrelation of BPSK signal (1-D)
% y = sigout;               % input signal (complex-valued)
% max_lag = 16;             % number of +/- time lags in correlation
% plotswitch = 1;           % generate plots 1->plots on, 0->plots off
% r_alpha1d = ralpha_1d(y,fcarr,rb,max_lag,fs,plotswitch);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Spectral Correlation  %%%%
%%%%   "Cyclic Spectrum"    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-Compute time-smoothed spectral correlation (2-D) "cyclic spectrum"
%-for selected baud rate
M = 32;                 % coherent averaging
L = 2;                 % incoherent averaging
alpha = rb;
plotswitch = 1;           % generate plots 1->plots on, 0->plots off
x = sigout(1:end);      % input signals (if x & y then compute cross-scd)
[~,~,~] = cyclospec_2d(fs,alpha,M,L,plotswitch,x);

% %-Compute time-smoothed spectral correlation (1-D) "cyclic spectrum" 
% %-for selected baud rate
% M = 32;                 % coherent averaging
% L = 16;                 % incoherent averaging
% alpha = rb;
% plotswitch = 1;           % generate plots 1->plots on, 0->plots off
% x = sigout(1:end);      % input signals (if x & y then compute cross-scd)
% [scd,~] = cyclospec_1d(fs,alpha,M,L,plotswitch,x);
