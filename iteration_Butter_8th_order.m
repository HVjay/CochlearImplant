file_name='different_filters/sound_quality_file4.wav'; % output signal of phase 1
final_filename='different_filters/sound_quality_file4_butter_8thOrder_lpf_kaiser.wav';
disp('ORDER 100')
[y,Fs]=audioread(file_name);
%sound(y,Fs);

[time_points,n]=size(y);
if (n==2)
    new_y=zeros(time_points,1);
    for k=1:time_points
        new_y(k,1)=y(k,1)+y(k,2);
    end
else
   new_y=y;
end

%sound(new_y,Fs);
    
   
[P,Q] = rat(Fs/16000);
y_resampled=resample(new_y,Q,P);


%create passband bank and set initial frequency range
passbank=zeros(8,length(y_resampled));
cutoff_1=100; % Start frequency%
cutoff_2=987.5; % End frequency%
k=size(passbank); 
rows=k(1);


for n=1:8
    N=8;
    Fs=16000;
    butter_iir=Butter_8th_order(N,Fs,cutoff_1,cutoff_2);
   
    tic
    filtered_signal=filter(butter_iir,y_resampled);
    toc
    
    %add filtered signal to passband bank and update frequency range
    passbank(n,:)=filtered_signal;
    cutoff_1=cutoff_2;
    cutoff_2=cutoff_2+987.5-1;
    
end

%create envelope and channels of each envelope
envelope=abs(passbank);
all_envelope=zeros(8,length(y_resampled));
k=size(all_envelope); 
rows=k(1);

for i=1:rows
    N=6;
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


[m,n]=size(y_resampled);
cosbank=zeros(8,length(y_resampled));


k=size(cosbank); 
rows2=k(1);
central_frequency=543.75;

for n=1:rows2
    [q,time_points]=size(all_envelope(n,:)); % time points from envelope signal
    t=time_points/Fs; % t = (Number of points total)/(number of points per second)
    T=1/Fs;
    num_of_points=linspace(0,t-T,time_points); % num_of_points = time points spread out from 0 to t

    %time = 0:T:t-T;
    %points_per_second=time_points/t;
    %central_frequency = (cutoff_1 + cutoff_2)/2;
    
    cos_function=cos(2*pi*central_frequency*num_of_points);
    %cos_function=cos_function.'; % Transpose matrix 
    
    % add filtered signal to cosband bank and update frequency range
    cosbank(n,:)= cos_function;
    central_frequency=central_frequency+937;
end

cutoff_1=100; % Start frequency%
cutoff_2=987.5; % End frequency%


modulated = cosbank.*all_envelope;
[m,n]=size(y_resampled);

final_signal=zeros(m,1);
for n=1:m
    final_signal(n)=sum(modulated(:,n));
end

max_val = max(abs(final_signal));
normalized = final_signal./max_val;


figure(4);
plot(normalized);
title('final signal');
xlabel('Time points');
ylabel('Gain (dB)');

figure(5)
plot(y_resampled)
title('original signal');
xlabel('Time points');
ylabel('Gain (dB)');

%filename='test_sounds/sound_quality_file2_400LPF.wav';
audiowrite(final_filename,normalized,Fs);