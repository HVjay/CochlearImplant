function Hd = kaiser_lpf(N,Fs,Fc,flag,Beta)
%KAISER_LPF Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.7 and Signal Processing Toolbox 8.3.
% Generated on: 03-Aug-2020 14:02:17

% FIR Window Lowpass filter designed using the FIR1 function.

% All frequency values are in Hz.
%Fs = 16000;  % Sampling Frequency

%N    = 8;        % Order
%Fc   = 400;      % Cutoff Frequency
%flag = 'scale';  % Sampling Flag
%Beta = 0.5;      % Window Parameter

% Create the window vector for the design algorithm.
win = kaiser(N+1, Beta);

% Calculate the coefficients using the FIR1 function.
b  = fir1(N, Fc/(Fs/2), 'low', win, flag);
Hd = dfilt.dffir(b);

% [EOF]
