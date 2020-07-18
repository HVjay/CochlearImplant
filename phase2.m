file_name='E_octave_new_new.wav'; % output signal of phase 1


[y, Fs] = audioread(file_name);
%{
[time_points,n]=size(y);
if (n==2)
    new_y=zeros(time_points,1);
    for k=1:time_points
        new_y(k,1)=y(k,1)+y(k,2);
    end
else
   new_y=y;
end

[P,Q] = rat(Fs/16000);
down_sampled=resample(new_y,Q,P);
filename='E_octave_new_new.wav';
audiowrite(filename,down_sampled,16000);
%}

passband=zeros(8,length(y));
Fc1=100;
Fc2=1000;
k=size(passband);
rows=k(1);
N=8;
flag='scale';
Beta=6;
fir_kaiser=FIR_kaiser(N,Fs,Fc1,Fc2,flag,Beta);
filtered_signal=filter(fir_kaiser,y);

    
L= length(filtered_signal);
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(filtered_signal,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
figure(2);
plot(f,2*abs(Y(1:NFFT/2+1)));
title('Single-Sided Amplitude Spectrum of y(t)');
xlabel('Frequency (Hz)');
ylabel('|Y(f)|');

%{
for n=1:rows
    N=8;
    flag='scale';
    Beta=6;
    fir_kaiser=FIR_kaiser(N,Fs,Fc1,Fc2,flag,Beta);
   % y_fft=fft(y);
    filtered_signal=filter(fir_kaiser,y);
    %plot(filtered_signal)
    
    L= length(filtered_signal);
    NFFT = 2^nextpow2(L); % Next power of 2 from length of y
    Y = fft(filtered_signal,NFFT)/L;
    plot(Y)
   % fft_new=fft(filtered_signal);
    %plot(fft_new);
end
%}