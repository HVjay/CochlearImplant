file_name='E_octave_new_new.wav'; % output signal of phase 1

[y, Fs] = audioread(file_name);

y_fft = fft(y);
L =length(y_fft);
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
figure(1);
plot(f,2*abs(Y(1:NFFT/2+1)));
title('Original y(t)');
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
filename='E_octave_new_new.wav';
audiowrite(filename,down_sampled,16000);
%}

%{
Fc1=100; % Start frequency%
Fc2=1000; % End frequency%
N=100;
flag='scale';
Beta=6;
Hd=FIR_Kaiser2(N,Fs,Fc1,Fc2,flag,Beta);
filtered_signal=filter(Hd,y);

L = length(filtered_signal);
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(filtered_signal,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
figure(2);
plot(f,2*abs(Y(1:NFFT/2+1)));
title('Filtered y(t)');
xlabel('Frequency (Hz)');
ylabel('|Y(f)|');
%}

passbank=zeros(8,length(y));
Fc1=100; % Start frequency%
Fc2=1000; % End frequency%
k=size(passbank); 
rows=k(1);


for n=1:rows
    N=100;
    flag='scale';
    Beta=6;
    
    
    fir_kaiser=FIR_Kaiser2(N,Fs,Fc1,Fc2,flag,Beta);
   % y_fft=fft(y);
   
    tic
    filtered_signal=filter(fir_kaiser,y);
    %plot(filtered_signal)
    toc
    
    passbank(n,:)=filtered_signal;
    Fc1=round(Fc2,-3);
    Fc2=round(Fc2+1000,-3)-1;
    
end


L = length(passbank(1, :));
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(passbank(1, :),NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
figure(3);
plot(f,2*abs(Y(1:NFFT/2+1)));
title(' Amplitude Spectrum of Filtered Sound waveform (100 - 1000)');
xlabel('Frequency (Hz)');
ylabel('Amplitude');


L = length(passbank(8, :));
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(passbank(8, :),NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
figure(4);
plot(f,2*abs(Y(1:NFFT/2+1)));
title(' Amplitude Spectrum of Filtered Sound waveform (7000 - 7999)');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

envelope=abs(passbank);

figure(10)
plot(envelope(1,:));

all_envelope=zeros(8,length(y));
for i=1:8
    Hd = lpf(Fs, 400, 10);
    all_envelope(i,:) = filter(Hd, envelope(i,:));
end

%[yupper,ylower]=envelope(

figure(5);
plot(all_envelope(1,:));

L = length(all_envelope(1, :));
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(all_envelope(1, :),NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
figure(6);
plot(f,2*abs(Y(1:NFFT/2+1)));
title(' Amplitude Spectrum of the envelope (100- 1000)');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
