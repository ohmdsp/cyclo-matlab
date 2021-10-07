function y = freq_shift(fnew,fold,x,fs);
% Author: D.R.Ohm
% Takes the input data sequence x and frequency shifts it by the 
% amount (fnew-fold)/fs. The shifted frequency sequence y is returned
%
% INPUTS:
% fnew - desired center frequency to shift data to
% fold - center frequency that data in x is currently at
% x - complex signal data sequence (1 X N) to frequency shift
% fs - sampling rate of x
%
% OUTPUTS:
% y - frequency shifted complex signal data stream
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% Calculate the difference in frequency between fnew and fold and normalize
% by fs
delta_f = (fnew-fold)./fs;
% Get the length of the signal stream x
N = length(x);
% Apply frequency shift
y = x.*exp(-1i*2*pi*delta_f*(0:(N-1)));