file_name='E_octave_new_new.wav'; % output signal of phase 1

[y, Fs] = audioread(file_name);

%create passband bank and set initial frequency range
passbank=zeros(8,length(y));
cutoff_1=100; % Start frequency%
cutoff_2=1000; % End frequency%
k=size(passbank); 
rows=k(1);


for n=1:rows
    N=25;
    flag='scale';
    Beta=6;
    
    fir_kaiser=FIR_Kaiser2(N,Fs,cutoff_1,cutoff_2,flag,Beta);
   
    %tic
    filtered_signal=filter(fir_kaiser,y);
    %toc
    
    %add filtered signal to passband bank and update frequency range
    passbank(n,:)=filtered_signal;
    cutoff_1=round(cutoff_2,-3);
    cutoff_2=round(cutoff_2+1000,-3)-1;
    
end

%plot of first and last channels
figure(1);
plot(passbank(1,:));
title('Filtered Sound waveform (100 - 1000 Hz)');
xlabel('Time points');
ylabel('Gain (dB)');

figure(2);
plot(passbank(8, :));
title('Filtered Sound waveform (7000 - 7999 Hz)');
xlabel('Time points');
ylabel('Gain (dB)');

%create envelope and channels of each envelope
envelope=abs(passbank);
all_envelope=zeros(8,length(y));
for i=1:8
    N=4;
    sampling_freq = 16000;
    cutoff_freq=400;
    Hd = butter_filter(N,sampling_freq,cutoff_freq);
    
    %add filtered envelope 
    all_envelope(i,:) = filter(Hd, envelope(i,:));
end

%plot envelopes
figure(3);
plot(all_envelope(1,:));
title('Envelope waveform (100 - 1000 Hz)');
xlabel('Time points');
ylabel('Gain (dB)');


figure(4);
plot(all_envelope(8, :));
title('Envelope waveform(7000 - 7999 Hz)');
xlabel('Time points');
ylabel('Gain (dB)');


