function Hd = least_squares_order2(N,Fs,Fpass1,Fstop1,Fpass2,Fstop2)
%LEAST_SQUARES_ORDER2 Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.7 and Signal Processing Toolbox 8.3.
% Generated on: 03-Aug-2020 00:35:06

% FIR least-squares Bandpass filter designed using the FIRLS function.

% All frequency values are in Hz.
%Fs = 16000;  % Sampling Frequency

%N      = 2;     % Order
%Fstop1 = 100;   % First Stopband Frequency
%Fpass1 = 100;   % First Passband Frequency
%Fpass2 = 1000;  % Second Passband Frequency
%Fstop2 = 1000;  % Second Stopband Frequency
Wstop1 = 1;     % First Stopband Weight
Wpass  = 1;     % Passband Weight
Wstop2 = 1;     % Second Stopband Weight

% Calculate the coefficients using the FIRLS function.
b  = firls(N, [0 Fstop1 Fpass1 Fpass2 Fstop2 Fs/2]/(Fs/2), [0 0 1 1 0 ...
           0], [Wstop1 Wpass Wstop2]);
Hd = dfilt.dffir(b);

% [EOF]
