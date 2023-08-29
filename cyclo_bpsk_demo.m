% Script: cyclo_demo.m
%
% Demo script showing common cycloprocessing functions and methods.
% 
% Author: drohm
%------------------------------------------------------------------------
%------------------------------------------------------------------------
clear all; close all;clc

%-Generate a simulated BPSK signal with AWGN
N  = 4*128;             % number of symbols
lmess = 4;              % time duration of signal (seconds)
fcarr = 100;              % carrier frequency (Hz)
samples = 8;            % samples per symbol
wgnvar = 0.2;          % variance of added noise (0 to 1)
plotswitch = 1;         % generate plots 1->plots on, 0->plots off
[sigout,fs,rb] = bpskgen(N,lmess,fcarr,samples,wgnvar,plotswitch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Cyclic Correlation    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-Compute non-conjugate cyclic autocorrelation of BPSK signal (2-D)
y = sigout;               % input signal (complex-valued)
max_lag = 16;             % number of +/- time lags in correlation
plotswitch = 1;           % generate plots 1->plots on, 0->plots off
r_alpha2d = ralpha_2d(y,fcarr,rb,max_lag,fs,plotswitch);

%-Compute cyclic autocorrelation of BPSK signal (1-D)
y = sigout;               % input signal (complex-valued)
max_lag = 16;             % number of +/- time lags in correlation
plotswitch = 1;           % generate plots 1->plots on, 0->plots off
r_alpha1d = ralpha_1d(y,fcarr,rb,max_lag,fs,plotswitch);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Spectral Correlation  %%%%
%%%%   "Cyclic Spectrum"    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-Compute time-smoothed spectral correlation (2-D) "cyclic spectrum"
%-for selected baud rate
M = 32;                 % coherent averaging
L = 16;                 % incoherent averaging
alpha = rb;
plotswitch = 1;           % generate plots 1->plots on, 0->plots off
x = sigout(1:end);      % input signals (if x & y then compute cross-scd)
[~,~,~] = cyclospec_2d(fs,alpha,M,L,plotswitch,x);

%-Compute time-smoothed spectral correlation (1-D) "cyclic spectrum" 
%-for selected baud rate
M = 32;                 % coherent averaging
L = 8;                 % incoherent averaging
alpha = rb;
plotswitch = 1;           % generate plots 1->plots on, 0->plots off
x = sigout(1:end);      % input signals (if x & y then compute cross-scd)
[scd,~] = cyclospec_1d(fs,alpha,M,L,plotswitch,x);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Spectral Coherence (TBD)   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %-Compute time-smoothed spectral coherence (2-D) "cyclic coherence"
% %-for selected baud rate
% M = 32;                 % coherent averaging
% L = 16;                 % incoherent averaging
% alpha = rb;
% x = sigout(1:end);      % input signals (if x & y then compute cross-scd)
% plotswitch = 1;
% [scd,~] = cyclo_coherence_2d(fs,alpha,M,L,plotswitch,x);


% %-Compute time-smoothed spectral coherence (1-D) ("cyclic spectrum" 
% %-for selected baud rate
% M = 32;                 % coherent averaging
% L = 16;                 % incoherent averaging
% alpha = rb;
% x = sigout(1:end);      % input signals (if x & y then compute cross-scd)
% [~,~] = cyclo_coherence_1d(fs,alpha,M,L,plotswitch,x);


