function [scd,freqx] = cyclodetect_fs(fs,alpha,M,plotswitch,x,y)
%
% Calculates the freqeuncy smoothed cyclic spectrum (also called the 
% Spectral Correlation Density (SCD) in liturature) for an input baud rate. 
% The cyclic correlation function and the spectral correlation function are
% Fourier transform pairs or cyclostationary signals, analogous to the 
% correlation and power spectraldensity pairs for stationary signals.
% An important feature of the second-order cyclosta-tionary statistics is 
% that they contain the phase information about the original signal, unlike
% the conventional PSD, which is a real-valued function.
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
t = 0:Ts:N*Ts-Ts;           % time index vector
Nfft = N;                   % FFT block size
window = hamming(N)';       % optional window function

for cnt=1:M
    upshft = exp(sqrt(-1)*pi*alpha*t);       
    dnshft = exp(-sqrt(-1)*pi*alpha*t);      
    u = x.*upshft;                    % shift signal up 1/2 cyclic freq
    v = y.*dnshft;                    % shift signal down 1/2 cyclic freq
    Syx = fftshift(fft(window.*u,Nfft).*conj(fft(window.*v,Nfft)));  
    Syx_tmp = [zeros(1,M/2), Syx ,zeros(1,M/2)];
      for ii = 1:N    
          Smyx(ii) = sum(Syx_tmp(ii:ii+M-1));   % moving average over M frames 
      end
    scdyx = Smyx;
    %cnt = cnt+1;
end
scd = 1/(M*N*Ts).*scdyx;       % scale the spectral correlation density

%-Generate Plots (optional)
if plotswitch == 1
    figure
    scd_mag = abs(scd);
    freqx = -fs/2:fs/(Nfft):fs/2-fs/(Nfft);     % frequency axis vector
    plot(freqx,scd_mag);
    grid on; 
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    title("Frequency-Smoothed Cyclic Spectrum at Baud Rate = ",+ alpha + " Hz" )     
end
