function r_alpha1d = ralpha_1d(y,fcarr,alpha,max_lag,fs,plotswitch)
% 
% Computes the 1-D cyclic autocorrelation function
%
% INPUT:
% y             - input signal (complex-valued)
% fcarr         - carrier freq of input signal
% alpha         - cyclic frequency
% lag           - number of +/- time lags in correlation
% fs            - sample rate of input signal
% plotswitch    - generate plots 1->plots on, 0->plots off
%
% OUTPUT:
% r_alpha       - non-conjugate cyclic correlation function
%
% Author: drohm
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
tau = -max_lag:max_lag;
 
yup = freq_shift(alpha/2,fcarr,y,fs);
ydn = freq_shift(-alpha/2,fcarr,y,fs);
r_alpha1d = correlation_sequence(max_lag,'biased',yup',ydn');
Cfmag_sq = abs(r_alpha1d);
    
if plotswitch == 1
    figure
    plot(tau,Cfmag_sq);
    axis tight;
    xlabel('tau (s)');ylabel('Magnitude');grid       
    title("1-D Cyclic Correlation at alpha = " + alpha + "Hz" )  
end