function [sigout,fs,rb] = bpskgen(N,lmess,fcarr,samples,wgnvar,plotswitch)
%
% Generate a BPSK signal in AWGN noise.
%
% INPUT:
% N             - number of symbols
% lmess         - time duration of message (seconds)
% fcarr         - carrier frequency (Hz)
% samples       - samples per symbol
% wgnvar        - variance of added noise (0 to 1)
% plotswitch    - generate plots 1->plots on, 0->plots off
%
% OUTPUT:
% sigout        - output bpsk signal
% fs            - sample frequency of bpsk output signal
% rb            - output symbol rate
%
% Author: drohm
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

T = lmess/N;            % symbol period (seconds)
Ts = T/samples;         % samples per symbol period
fs = 1/Ts;              % sample rate (samples/sec)
rb = N/lmess;           % symbol rate of BPSK output signal
td=(0:Ts:lmess-Ts)';    % Time vector
%tds=(0:Ts/samples:lmess-Ts/samples)';    % Time vector

fn = fs/2;              % Nyquist frequency
if fcarr > fn/2
    disp(' ')
    disp('Carrier frequency must be lower than 1/4 Nyquist for this function')
    disp(['Nyquist freuency is ' , num2str(fn)]) 
    error ('Carrier frequency too high!')
end

%-Generate BPSK symbols
data = sign(randn(N,1))';               % generate N random binary symbols, +/- 1(notice transpose')
symbols = ones(T/Ts,1)*data;            % upsampled symbols 
%symbols = ones(N,1)*data;               % non-upsampled symbols
Isymbols = symbols(:)';                 % form Inphase (I) data stream
Qsymbols = zeros(1,length(Isymbols));   % form Quadrature (Q) data stream (BPSK=0)
Csymbols = Isymbols + 1i*Qsymbols;       % complex baseband signal (no noise)

xsig = Csymbols;
% %-Design square-root-raised-cosine (SRRC) filter. 
% rolloff = .25;
% span = 6; 
% b = rcosdesign(rolloff, span, samples, 'normal');
% %fvtool(h, 'Analysis', 'impulse')   % Visualize the filter
% 
% %-Filter the data at the transmitter
% x = upfirdn(Csymbols, b, 1);
% xsig = x(1:end-samples*span);

%-Generate AWGN noise
sigma_w = wgnvar;                       % noise variance from input
Iwn = sigma_w*randn(1,length(xsig));
Qwn = sigma_w*randn(1,length(xsig));
wn = Iwn + 1i*Qwn;                % complex white gaussian noise

%-Add AWGN to complex baseband BPSK signal
bpsk_s = xsig + wn;

%-Generate carrier wave
twopi_fc_t = td*2*pi*fcarr;     % argument for carrier sinusoid 
a = 1;                          % amplitude
phi = 0;                        % initial phase
cs_t = a*exp(-1i*twopi_fc_t + phi);
%cs_t = a*cos(twopi_fc_t + phi); % carrier signal
cs_t=cs_t';   % transpose to column vector

%-Mix I/Q data stream with carrier
bpsk_c = cs_t.*bpsk_s;          % mix symbols with carrier
bpsk_c = .7*bpsk_c;             % reduce gain to keep output at +/- one

%-Final BPSK output signal
sigout = bpsk_c;

%-Generate Plots (optional)
if plotswitch == 1
    figure
    subplot(4,1,1)
    plot(td,Isymbols)
    ylim([-2 2]); grid
    xlabel('Time')
    ylabel('Amplitude')
    title('I-Symbol Data')

    subplot(4,1,2)
    plot(td,Qsymbols)
    ylim([-2 2]); grid
    xlabel('Time')
    ylabel('Amplitude')
    title('Q-Symbol Data')

    subplot(4,1,3)
    plot(td,real(bpsk_s))
    ylim([-2 2]); grid
    xlabel('Time')
    ylabel('Amplitude')
    title('BPSK Waveform with AWGN')
   
    %-Take FFT of BPSK modulated carrier
    y = bpsk_c;
    nfft = 2.^(ceil(log(length(y))/log(2)));
    ffty = fft(y,nfft);                         % pad with zeros and take fft
    ffty_mag = abs(fftshift(ffty));
    freq_limit = fs/2;
    freq= -freq_limit:fs/nfft:freq_limit-fs/nfft;
    
    %-Plot frequency domain
    subplot(4,1,4)
    plot(freq,10*log10(ffty_mag.^2))
    xlim([-freq_limit,freq_limit]);ylim([-10 80]);grid
    xlabel('FREQUENCY(Hz)');ylabel('DB');
    title('BPSK Signal Spectrum')
    
%     %-Plot BPSK symbol constellation
%     scatterplot(bpsk_s); grid       % plot bpsk symbol constellation
%     title('BPSK Signal Constellation')
    
end  % end plots section
