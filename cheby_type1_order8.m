function Hd = cheby_type1_order8(N,Fs,Fpass1,Fpass2)
%CHEBY_TYPE1_ORDER8 Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.7 and Signal Processing Toolbox 8.3.
% Generated on: 02-Aug-2020 14:56:45

% Chebyshev Type I Bandpass filter designed using FDESIGN.BANDPASS.

% All frequency values are in Hz.
%Fs = 16000;  % Sampling Frequency

%N      = 8;     % Order
%Fpass1 = 100;   % First Passband Frequency
%Fpass2 = 1000;  % Second Passband Frequency
Apass  = 1;     % Passband Ripple (dB)

% Construct an FDESIGN object and call its CHEBY1 method.
h  = fdesign.bandpass('N,Fp1,Fp2,Ap', N, Fpass1, Fpass2, Apass, Fs);
Hd = design(h, 'cheby1');

% [EOF]
