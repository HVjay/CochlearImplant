file_name='E_octave_new_new.wav'; % output signal of phase 1


[y, Fs] = audioread(file_name);
figure (3);
plot (y);

y_fft = fft(y);
L =length(y_fft);
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
figure(1);
plot(f,2*abs(Y(1:NFFT/2+1)));
title('Single-Sided Amplitude Spectrum of original y(t)');
xlabel('Frequency (Hz)');
ylabel('|Y(f)|');

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
filename='running_new_new.wav';
audiowrite(filename,down_sampled,16000);
%}

passbank=zeros(8,length(y));
Fc1=100; % Start frequency%
Fc2=1000; % End frequency%
k=size(passbank); 
rows=k(1);
N=100;
flag='scale';
Beta=6;
Hd=FIR_kaiser(N,Fs,Fc1,Fc2,flag,Beta);
filtered_signal=filter(Hd,y);

L = length(filtered_signal);
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(filtered_signal,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
figure(2);
plot(f,2*abs(Y(1:NFFT/2+1)));
title('Single-Sided Amplitude Spectrum of filtered y(t)');
xlabel('Frequency (Hz)');
ylabel('|Y(f)|');


for n=1:rows
    N=100;
    flag='scale';
    Beta=6;
    disp(Fc1);
    fir_kaiser=FIR_kaiser(N,Fs,Fc1,Fc2,flag,Beta);
   % y_fft=fft(y);
    filtered_signal=filter(fir_kaiser,y);
    passbank(n,:)=filtered_signal;
    Fc1=Fc2;
    Fc2=Fc2+1000-1;
    %plot(filtered_signal)
  
   % fft_new=fft(filtered_signal);
    %plot(fft_new);
end


envelope=abs(passbank);
all_envelope=zeros(8,length(y));
for i=1:8
    Hd = lpf(Fs, 400, 10);
    all_envelope(i,:) = filter(Hd, envelope(i,:));
end

L = length(filtered_signal);
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(filtered_signal,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
figure(2);
plot(f,2*abs(Y(1:NFFT/2+1)));
title('Single-Sided Amplitude Spectrum of filtered y(t)');
xlabel('Frequency (Hz)');
ylabel('|Y(f)|');
plot(all_envelope(1, :)), xlabel('Sample Number'),ylabel('Amplitude'),title('Filtered Sound waveform (100 - 460)');
plot(all_envelope(22, :)), xlabel('Sample Number'),ylabel('Amplitude'),title('Filtered Sound waveform (7640 -8000)');



