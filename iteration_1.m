file_name='2.3_ambulance.wav'; % output signal of phase 1

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
passbank=zeros(22,length(y_resampled));
cutoff_1=100; % Start frequency%
cutoff_2=459.1; % End frequency%
k=size(passbank); 
rows=k(1);


for n=1:22
    N=100;
    flag='scale';
    Beta=6;
    
    fir_kaiser=FIR_Kaiser2(N,Fs,cutoff_1,cutoff_2,flag,Beta);
   
    %tic
    filtered_signal=filter(fir_kaiser,y_resampled);
    %toc
    
    %add filtered signal to passband bank and update frequency range
    passbank(n,:)=filtered_signal;
    disp(cutoff_1)
    disp(cutoff_2)
    cutoff_1=cutoff_2;
    cutoff_2=cutoff_2+359.1-1;
    
end

%create envelope and channels of each envelope
envelope=abs(passbank);
all_envelope=zeros(22,length(y_resampled));
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
cosbank=zeros(22,length(y_resampled));


k=size(cosbank); 
rows2=k(1);
central_frequency=279.6;

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
    central_frequency=central_frequency+359.1;
end

cutoff_1=100; % Start frequency%
cutoff_2=500; % End frequency%


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

filename='2.3_ambulance_new.wav';
audiowrite(filename,normalized,Fs);