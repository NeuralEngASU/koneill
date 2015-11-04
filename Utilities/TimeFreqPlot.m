% x = testData.Data(20,125000:350000);
% x = tmpdata(:,1);
% x = tmpData1;
% x = detrend(data(1,1:300000));%data(:,1,1);
x = double(data(:,1));
figure;
% Fs = testData.MetaTags.SamplingFreq; % sampling frequency
Fs = 500;
L = length(x); % length of signal
t = 0:1/Fs:(L-1)/Fs; % time base
NFFT = 2^nextpow2(L);
f = Fs/2*linspace(0,1,NFFT/2+1); % single sided spectrum

% high-pass filter 
order = 3;
Fc = 10; % cutoff frequency
[z,p,k] = butter(order,Fc/(Fs/2),'high');
[SOS,G] = zp2sos(z,p,k);% convert to SOS structure to use filter analysis tool

x_filt = filtfilt(SOS,G,x);

tic;
[pxx,f] = pmtm(x,9,NFFT,Fs);
toc

figure;
ax1 = subplot(3,1,1); % filtered time-voltage.
plot(t,x,'k');% time-voltage.
title('Broadband')
xlabel('Time (seconds)')
ylabel('Volts')
ax2 = subplot(3,1,2);
plot(t,x_filt,'k');
title('High-Pass Filtered');
xlabel('Time (seconds)');
ylabel('Volts');
subplot(3,1,3);
loglog(f,pxx);
title('Multi-taper Spectrum');
xlabel('Frequency (Hz)')
ylabel('Power')
linkaxes([ax1,ax2],'x');


% multi-taper power spectrum with 95% confidence bounds.
[pxx,f,pxxc] = pmtm(x_filt,9,NFFT,Fs,'ConfidenceLevel',0.95);
title('Multitaper PSD Estimate with 95%-Confidence Bounds')