function Hd = FIR_kaiser(N,Fs,Fc1,Fc2,flag,Beta)
%FIR_KAISER Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.7 and Signal Processing Toolbox 8.3.
% Generated on: 18-Jul-2020 17:56:04

% FIR Window Bandstop filter designed using the FIR1 function.

% All frequency values are in Hz.
%Fs = 16000;  % Sampling Frequency

%{
N    = 8;        % Order
Fc1  = 100;      % First Cutoff Frequency
Fc2  = 200;      % Second Cutoff Frequency
flag = 'scale';  % Sampling Flag
Beta = 6;        % Window Parameter
%}

% Create the window vector for the design algorithm.
win = kaiser(N+1, Beta);

% Calculate the coefficients using the FIR1 function.
b  = fir1(N, [Fc1 Fc2]/(Fs/2), 'stop', win, flag);
Hd = dfilt.dffir(b);

% [EOF]