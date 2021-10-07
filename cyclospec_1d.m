function [scd,freqx] = cyclospec_1d(fs,alpha,M,L,plotswitch,x,y)
%
% Calculates the (1D) time-smoothed non-conjugate cyclic spectrum 
% (Spectral Correlation) for an input baud rate. 
% The cyclic correlation function and the spectral correlation function are
% Fourier transform pairs for cyclostationary signals, analogous to the 
% correlation and power spectraldensity pairs for stationary signals.
% An important feature of the second-order cyclostationary statistics is 
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
% L             -   Number of incoherent averaging blocks
% plotswitch    -   generate plots 1->plots on, 0->plots off
% x and y       -   Input signals (if both then compute cross-spectra)
% 
% OUTPUTS:
% scd           -   Time-smoothed spectral correlation (non-conjugate)
% freqx         -   Indices for frequency axis
%
% Author: drohm
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

if nargin < 7
    y = x;
end
N = length(x);                  % number of samples for input signal
Ts = 1/fs;                      % sample period
Nfft = 4*1024;                  % FFT block size
recordsize = floor(N/L);        % signal record size for incoherent averaging
blocksize = floor(recordsize/M);% block size given smoothing factor M for coherent averaging
tb = 0:Ts:blocksize*Ts-Ts;      % time index vector for block size

window = 1;%hamming(blocksize)';            % optional window function
scd = zeros(1,Nfft);
for j = 1:L   % incoherent averaging
    xr = x((j-1)*recordsize+1:j*recordsize);
    yr = y((j-1)*recordsize+1:j*recordsize);

    scdyx = zeros(1,Nfft);
    for i = 1:M     % coherent integration
        upshft = exp(1i*pi*alpha*tb);       
        dnshft = exp(-1i*pi*alpha*tb);      
        u = xr((i-1)*blocksize+1:i*blocksize).*upshft;   % shift signal up 1/2 cyclic freq
        v = yr((i-1)*blocksize+1:i*blocksize).*dnshft;   % shift signal down 1/2 cyclic freq
        scdyx = scdyx + fftshift(fft(window.*u,Nfft).*conj(fft(window.*v,Nfft)));    
    end         
    scd_tmp = scdyx;            % spectral correlation   
end
scd = scd + abs(scd_tmp);

freqx = -fs/2:fs/(Nfft):fs/2-fs/(Nfft);     % frequency axis vector

%-Generate Plots (optional)
if plotswitch == 1
    figure
    subplot(2,1,1)
    plot(freqx,scd);
    grid on;axis tight
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    title("Cyclic Spectrum (Spectral Correlation) at alpha = " + alpha + "Hz" ) 
    subplot(2,1,2)
    plot(freqx,180/pi.*angle(scd_tmp), 'r')
    grid on;axis tight;ylim([-180,180])
    xlabel('Frequency (Hz)')
    ylabel('Phase Angle (deg)')
    title("Cyclic Phase Spectrum at alpha = " + alpha + "Hz" ) 


end
