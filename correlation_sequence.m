function r=correlation_sequence(max_lag,bias,x,y)

% Computes the unbiased/biased auto correlation sequence (ACS) or cross
% correlation sequence (CCS) estimate. Most computationally efficient if
% data vector length is a power of 2.  Returns length 2*max_lag + 1
% correlation sequence in a column vector.  Zero lag of correlation falls
% in middle of the sequence at element (row) max_lag + 1.  'biased'
% estimates scaled by 1/length(x); 'unbiased' estimates scaled by
% 1/[length(x) - lag index].
%
%     ACS:  r=correlation_sequence(max_lag,bias,x)
%     CCS:  r=correlation_sequence(max_lag,bias,x,y)
%
% max_lag -- maximum time index for estimated ACS or CCS
% bias    -- select: 'unbiased' or 'biased' estimates
% x       -- column vector of data samples from channel x
% y       -- column vector of data samples from channel y (CCS only)
% r       -- ACS or CCS vector ordered from -max_lag to +max_lag;
%               zero lag is element max_lag + 1

nx = length(x);
if max_lag > nx, error('Parameter MAX_LAG > length(x).'), end
if nargin < 4
    temp = fft([x;zeros(size(x))]);
    r = fftshift(ifft(real(temp).^2 + imag(temp).^2));
else
    ny = length(y);
    if ny ~= nx, error('Data vectors x and y not of equal length.'), end
    r = fftshift(ifft(fft([x;zeros(size(x))]).*conj(fft([y;zeros(size(y))]))));  
end
r = r(nx+1-max_lag:nx+1+max_lag);               % truncate to requested max_lag
if strcmp(bias,'biased')
    r = r/nx;                                   % eqs (5.13),(5.19)
else                                            % default is 'unbiased'
    r = r./(nx-abs(max_lag-(0:2*max_lag)))';    % eq (5.9)
end
%
