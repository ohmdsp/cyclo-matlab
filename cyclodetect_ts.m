function [scd,freqx] = cyclodetect_ts(fs,alpha,M,plotswitch,x,y)
%
% Calculates the time-smoothed cyclic spectrum (also called the 
% Spectral Correlation Density (SCD) in liturature) for an input baud rate. 
% The cyclic correlation function and the spectral correlation function are
% Fourier transform pairs or cyclostationary signals, analogous to the 
% correlation and power spectraldensity pairs for stationary signals.
% An important feature of the second-order cyclostationary statistics is 
% that they contain the phase information about the original signal, unlike
% the conventional PSD, which is a real-valued function.
% Note - important to average the magnitudes, not the complex spectrums
% (since cannot assume to be coherent unless phase is maintained accross
% block sizes).
%
% Notes: 
% Integration time T should not be greater than coherence (reciprocal
% of bandwidth). 
% M should be choosen so that deltaf*T >> 1 (not to exceed bandwidth).
%
% INPUTS:
% fs            -   Sample frequency for input signals
% alpha         -   baud rate for analysis
% M             -   Frequency smoothing factor
% plotswitch    -   generate plots 1->plots on, 0->plots off
% x and y       -   Input signals (if both then compute cross-spectra)
% 
% OUTPUTS:
% scd           -   Frequency smoothed spectral correlation density
% freqx         -   Indices for frequency axis
%
% Author: drohm
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

if nargin < 6
    y = x;
end
N = length(x);              % number of samples for input signal
Ts = 1/fs;                  % sample period
Nfft = N;                   % FFT block size
blocksize = floor(N/M);    % block size given smoothing factor M
t = 0:Ts:blocksize*Ts-Ts;           % time index vector
 
window = hamming(blocksize)';            % optional window function
scd = zeros(1,Nfft);

for i = 1:M
    upshft = exp(sqrt(-1)*pi*alpha*t);       
    dnshft = exp(-sqrt(-1)*pi*alpha*t);      
    u = x((i-1)*blocksize+1:i*blocksize).*upshft;   % shift signal up 1/2 cyclic freq
    v = y((i-1)*blocksize+1:i*blocksize).*dnshft;   % shift signal down 1/2 cyclic freq
    scdyx = fftshift(fft(window.*u,Nfft).*conj(fft(window.*v,Nfft)));
    scd = scd + 1/(M*N*Ts).*abs(scdyx);            % scale the spectral correlation density
end         
freqx = -fs/2:fs/(Nfft):fs/2-fs/(Nfft);     % frequency axis vector

%-Generate Plots (optional)
if plotswitch == 1
    figure
    scd_mag = scd;
    plot(freqx,scd_mag);
    grid on; 
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    title("Time-Smoothed Cyclic Spectrum at Baud Rate = ",+ alpha + " Hz" )     
end