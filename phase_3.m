file_name='restuarant_new.wav'; % output signal of phase 1

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

%-------------------------------------------------------------------------%
% PHASE 3 STARTS HERE
%-------------------------------------------------------------------------%
%create cosine function bank and set initial frequency range

cosbank=zeros(8,length(y));
cutoff_1=100; % Start frequency%
cutoff_2=1000; % End frequency%
k=size(cosbank); 
rows2=k(1);

for n=1:rows2
    [q,time_points]=size(all_envelope(n,:)); % time points from envelope signal
    t=time_points/Fs; % t = (Number of points total)/(number of points per second)

    num_of_points=linspace(0,t,time_points); % num_of_points = time points spread out from 0 to t

    points_per_second=time_points/t;
    central_frequency = (cutoff_1 + cutoff_2)/2;
    
    cos_function=cos(2*pi*central_frequency*num_of_points);
    %cos_function=cos_function.'; % Transpose matrix 
    
    % add filtered signal to cosband bank and update frequency range
    cosbank(n,:)= cos_function;
    cutoff_1=round(cutoff_2,-3);
    cutoff_2=round(cutoff_2+1000,-3)-1;

end

%{
points_per_5ms= floor((points_per_second/1000)*5);
figure(1);
plot(cosbank(1, 1:points_per_5ms));
title('Cosbank channel 1');
xlabel('Time points');
ylabel('Gain (dB)');
%}

modulated=zeros(8,length(y));
k=size(modulated); 
rows2=k(1);

cutoff_1=100; % Start frequency%
cutoff_2=1000; % End frequency%
for n=1:rows2
    %cosbank_channel_transposed=cosbank(n,:).';
    %modulated_signal=cosbank_channel_transposed* all_envelope(n,:);
    %negative_channel=-(all_envelope(n,:));
    modulated_signal_positive=cosbank(n,:).*all_envelope(n,:);
    
    %central_frequency = (cutoff_1 + cutoff_2)/2;
    %modulated_signal=ammod(all_envelope(n,:),central_frequency,16000);
    
    %modulated_channel_negative=cosbank(n,:).*negative_channel;
    
    modulated(n,:)=modulated_signal_positive;
    cutoff_1=round(cutoff_2,-3);
    cutoff_2=round(cutoff_2+1000,-3)-1;
end    

%{
figure(2);
plot(all_envelope(2, :));
title('rectified channel 2');
xlabel('Time points');
ylabel('Gain (dB)');

figure(3);
plot(modulated(2, :));
title('modulated channel 2');
xlabel('Time points');
ylabel('Gain (dB)');
%}

final_signal=zeros(1,length(y));
for n=1:rows2
    final_signal=final_signal+modulated(n,:);
end

final_abs=abs(final_signal);
final_abs_max=max(final_abs);
final_signal=final_signal./final_abs_max;


figure(4);
plot(final_signal);
title('final signal');
xlabel('Time points');
ylabel('Gain (dB)');

figure(5)
plot(y)
title('original signal');
xlabel('Time points');
ylabel('Gain (dB)');

filename='restuarant_phase3.wav';
audiowrite(filename,final_signal,Fs);