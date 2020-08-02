function Hd = Butter_8th_order(N,Fs,Fc1,Fc2)
%BUTTER_8TH_ORDER Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.7 and Signal Processing Toolbox 8.3.
% Generated on: 02-Aug-2020 14:35:25

% Butterworth Bandpass filter designed using FDESIGN.BANDPASS.

% All frequency values are in Hz.
%Fs = 48000;  % Sampling Frequency

%N   = 8;     % Order
%Fc1 = 100;   % First Cutoff Frequency
%Fc2 = 1000;  % Second Cutoff Frequency

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.bandpass('N,F3dB1,F3dB2', N, Fc1, Fc2, Fs);
Hd = design(h, 'butter');

% [EOF]
