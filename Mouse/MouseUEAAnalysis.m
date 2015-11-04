%% 

% Reads in data file.
[t,amps,data,aux] = read_intan_data;


span = 2000:12000;

pk2pk = zeros(1,16);
avgDC = zeros(1,16);
for k = 1:16
    pk2pk(k) = max(data(span,k)) - min(data(span, k));
    avgDC(k) = mean(data(span,k));
end

figure (1)

[AX,H1,H2] = plotyy(1:16, pk2pk./1000, 1:16, avgDC./1000);
grid on
title('Mouse Optogenetics: No UEA Connected')
xlabel('Channel')
% ylabel(sprintf('Peak-to-Peak Voltage (mV)'))
set(get(AX(1),'Ylabel'),'String','Peak-to-Peak Voltage (mV)') 
set(get(AX(2),'Ylabel'),'String','DC Offset Voltage (mV)') 
set(AX(1),'XLim', [1, 16]) 
set(AX(2),'XLim', [1, 16]) 


%%

% Reads in data file.
[t,amps,data,aux] = read_intan_data;


span = 2000:12000;

pk2pk = zeros(1,16);
avgDC = zeros(1,16);
for k = 1:16
    pk2pk(k) = max(data(span,k)) - min(data(span, k));
    avgDC(k) = mean(data(span,k));
end

figure (2)

[AX,H1,H2] = plotyy(1:16, pk2pk./1000, 1:16, avgDC./1000);
grid on
title('Mouse Optogenetics: UEA Connected, No Saline')
xlabel('Channel')
% ylabel(sprintf('Peak-to-Peak Voltage (mV)'))
set(get(AX(1),'Ylabel'),'String','Peak-to-Peak Voltage (mV)') 
set(get(AX(2),'Ylabel'),'String','DC Offset Voltage (mV)') 
% set(get(AX(1),'XLim'),[1, 16]) 
% set(get(AX(2),'XLim'), [1, 16]) 

set(AX(1),'XLim', [1, 16]) 
set(AX(2),'XLim', [1, 16]) 

%%
% Reads in data file.
[t,amps,data,aux] = read_intan_data;


span = 2000:12000;

pk2pk = zeros(1,16);
avgDC = zeros(1,16);
for k = 1:16
    pk2pk(k) = max(data(span,k)) - min(data(span, k));
    avgDC(k) = mean(data(span,k));
end

figure (2)

[AX,H1,H2] = plotyy(1:16, pk2pk./1000, 1:16, avgDC./1000);
grid on
title('Mouse Optogenetics: UEA in Saline, Ground and Referenece in Saline')
xlabel('Channel')
% ylabel(sprintf('Peak-to-Peak Voltage (mV)'))
set(get(AX(1),'Ylabel'),'String','Peak-to-Peak Voltage (mV)') 
set(get(AX(2),'Ylabel'),'String','DC Offset Voltage (mV)') 
% set(get(AX(1),'XLim'),[1, 16]) 
% set(get(AX(2),'XLim'), [1, 16]) 

set(AX(1),'XLim', [1, 16]) 
set(AX(2),'XLim', [1, 16]) 

%%

channel = 4;

[t,amps,data,aux] = read_intan_data;

figure (channel)

Fs = 25000;
span = (1:length(data(:,channel)))';
% span = length(data(:,1));

NFFT = 2^nextpow2(length(data(:,channel))); % Next power of 2 from length of y

f = Fs/2*linspace(0,1,NFFT/2+1);

Y = zeros(NFFT , 1);
Yavg = zeros(NFFT , 1);
for k = 1:16
%     pk2pk(k) = max(data(span,k)) - min(data(span, k));
%     avgDC(k) = mean(data(span,k));

    Y(:,k) = fft(data(span, k),NFFT)/length(data(:,channel));
    
    Yavg = [Yavg + Y(:,k)];

end

Yavg = Yavg/16;

% Plot single-sided amplitude spectrum.
% plot(f(1:500),2*abs(Yavg(1:500)))
% plot(f(1:500),2*abs(Y(1:500,4)))

plot(f,2*abs(Y(1:length(f),channel)))
title(sprintf('Single-Sided Amplitude Spectrum of Y(t) [Channel: %d]\n with input of 1 mV at 1500 Hz', channel))
% title(sprintf('Single-Sided Amplitude Spectrum of Y(t) [Channel: %d]\n UEA in Saline, No Input', channel))
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
xlim([0, 2000])


% plot(t(span), data(span, 1))
