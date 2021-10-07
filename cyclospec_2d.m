function [scd,freqx,alphax] = cyclospec_2d(fs,alpha,M,L,plotswitch,x,y)
%
% Calculates the (2D) time-smoothed cyclic spectrum (Spectral Correlation)
% for an input baud rate. 
% The cyclic correlation function and the spectral correlation function are
% Fourier transform pairs for cyclostationary signals, analogous to the 
% correlation and power spectraldensity pairs for stationary signals.
% An important feature of the second-order cyclostationary statistics is 
% that they contain the phase information about the original signal, unlike
% the conventional PSD, which is a real-valued function.
%
% Note: alpha/2 shifts should be multiples of 1/M. If not, zero pad until
% it is, or very close. 
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
Nfft = 2*1024;                  % FFT block size
recordsize = floor(N/L);        % signal record size for incoherent averaging
blocksize = floor(recordsize/M);% block size given smoothing factor M for coherent averaging
tb = 0:Ts:blocksize*Ts-Ts;      % time index vector for block size

window = 1;%hamming(blocksize)';            % optional window function
stepsz = alpha;                           % evaluation step for plots
range = 4;
%alphax = -range*alpha:stepsz:range*alpha;
alphax = 0:stepsz:range*alpha;
%alphax = -range*alpha:stepsz:0;
freqx = -fs/2:fs/(Nfft):fs/2-fs/(Nfft);     % frequency axis vector
scd = zeros(length(alphax),length(freqx));
index = 1;
%i = range*alpha:-1:0
for nn = range*alpha:-1*stepsz:0
%for nn = 0:stepsz:range*alpha
    alpha_shift = nn-range*alpha;
    
    for j = 1:L   % incoherent averaging
        xr = x((j-1)*recordsize+1:j*recordsize);
        yr = y((j-1)*recordsize+1:j*recordsize);

        scdyx = zeros(1,Nfft);
        for i = 1:M     % coherent integration
            upshft = exp(1i*pi*alpha_shift*tb);       
            dnshft = exp(-1i*pi*alpha_shift*tb);      
            u = xr((i-1)*blocksize+1:i*blocksize).*upshft;   % shift signal up 1/2 cyclic freq
            v = yr((i-1)*blocksize+1:i*blocksize).*dnshft;   % shift signal down 1/2 cyclic freq
            scdyx = scdyx + fftshift(fft(window.*u,Nfft).*conj(fft(window.*v,Nfft)));
        end         
        scd_tmp = 1/M.*scdyx;                    % spectral correlation   
    end
    scd(index,:) = 1/L.*abs(scd_tmp);
    index = index + 1;
end


if plotswitch == 1
    figure
    h = waterfall(freqx,alphax,scd);
    axis tight;view([-150 20])
    set(h, 'FaceColor', '[ 0.9482 0.8157 0.6443]');   
    set(h, 'EdgeColor', 'k');
    h.FaceAlpha = 0.90;  
    xlabel('freq (Hz)');ylabel('alpha (Hz)');grid on      
    title("2-D Cyclic Spectrum (Spectral Correlation)" )      
end


