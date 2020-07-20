function Hd = lpf(Fs, Fc, N)
%LPF Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.8 and Signal Processing Toolbox 8.4.
% Generated on: 14-Jul-2020 18:44:15

% Butterworth Lowpass filter designed using FDESIGN.LOWPASS.

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.lowpass('N,F3dB', N, Fc, Fs);
Hd = design(h, 'butter');

% [EOF]
end