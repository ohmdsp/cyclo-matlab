function r_alpha2d = ralpha_2d(y,fcarr,alpha,max_lag,fs,plotswitch)
% 
% Computes the 2-D non-conjugate cyclic autocorrelation function
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
% r_alpha       - cyclic correlation function
%
% Author: drohm
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
tau = -max_lag:max_lag;
range = 4;
%alphax = -range*alpha:range*alpha;
alphax = 0:range*alpha;

index = 1;
for i = range*alpha:-1:0
    alpha_shift = i-range*alpha;
    yup = freq_shift(alpha_shift/2,fcarr,y,fs);
    ydn = freq_shift(-alpha_shift/2,fcarr,y,fs);
    r_alpha2d = correlation_sequence(max_lag,'biased',yup',ydn');
    Cfmag_sq(index,:) = abs(r_alpha2d);
    index = index + 1;   
end

if plotswitch == 1
    figure
    stepsz = alpha;
    h = waterfall(tau,alphax(1:stepsz:end), Cfmag_sq(1:stepsz:end,:));
    axis tight;view([-150 20])
    set(h, 'FaceColor', '[ 0.6443 0.8157 0.9482]');   
    set(h, 'EdgeColor', 'k');
    h.FaceAlpha = 0.90;  
    xlabel('tau (s)');ylabel('alpha (Hz)');grid on       
    title("2-D Cyclic Correlation" )     
end