%% Covariance
load('C:\Users\admin\Desktop\Cheetah-acquisition files\Raw data files\2013-05-07_08-40-59\exp2013-05-07_08-40-59_CSC.mat')

% Grabs the sample data

for i=1:4
    test(:,i) = CSCData.samples{i}';
end % END FOR

test = test./10e2;

C = cov(test);

%% Correlation Coeficient
% A normalized magnitude version of the correlation matrix. Where 1.0 =
% perfect correlation, and 0 means no correlation


load('C:\Users\admin\Desktop\Cheetah-acquisition files\Raw data files\2013-05-07_08-40-59\exp2013-05-07_08-40-59_CSC.mat')

% Grabs the sample data

for i=1:4
    test(:,i) = CSCData.samples{i}';
end % END FOR

% test = test./10e2; % Reduces the magnitude

R = corrcoef(test);
% Rij = Cij/sqrt(Cii*Cjj)

%% Plot CorrCoef
% Assume that the 4 electrodes are physically alligned like:
%  1 2
%  3 4

numSquare = 4*4; % Each electrode will have four plots associated with it.

colormap(jet); % Sets the current color map to jet

mapColor = jet; % Contains the RGB valuse for the jet color map.

sqPchCoord = [4,0; 4,2; 2,0; 2,2]; % Coordinates for the large patches

hold on
for i=1:4

    for j=1:4
        
        squarePatchxx = [sqPchCoord(i,2), sqPchCoord(i,2)+1, sqPchCoord(i,2)+1, sqPchCoord(i,2)];
        squarePatchyy = [sqPchCoord(i,1), sqPchCoord(i,1), sqPchCoord(i,1)-1, sqPchCoord(i,1)-1];
        
        switch j
            case 1
                squarePatchxx = squarePatchxx;
                squarePatchyy = squarePatchyy;
            case 2
                squarePatchxx = squarePatchxx + 1;
                squarePatchyy = squarePatchyy;
            case 3
                squarePatchxx = squarePatchxx;
                squarePatchyy = squarePatchyy - 1;
            case 4
                squarePatchxx = squarePatchxx + 1;
                squarePatchyy = squarePatchyy - 1;
            otherwise
        end % END SWITCH
            
        smallPatch = patch(squarePatchxx, squarePatchyy, 1);
        
        set(smallPatch, 'FaceColor', mapColor(round(R(i,j)/R(i,i) * 64),:))
        set(smallPatch, 'EdgeColor', 'none')
        
    end % END FOR
end % END FOR

colorbar('EastOutside')
colorbar('YTickLabel',...
    linspace(0,1,11))

for i=1:4
    squarePatchxx = [sqPchCoord(i,2), sqPchCoord(i,2)+2, sqPchCoord(i,2)+2, sqPchCoord(i,2)];
    squarePatchyy = [sqPchCoord(i,1), sqPchCoord(i,1), sqPchCoord(i,1)-2, sqPchCoord(i,1)-2];
    
    largePatch = patch(squarePatchxx, squarePatchyy, 1);
    
    set(largePatch, 'FaceColor', 'none')
    set(largePatch, 'EdgeColor', 'k')
end % END FOR


hold off

%% Negative correlation test

load('C:\Users\admin\Desktop\Cheetah-acquisition files\Raw data files\2013-05-07_08-40-59\exp2013-05-07_08-40-59_CSC.mat')

% Grabs the sample data

for i=1:4
    test(:,i) = CSCData.samples{i}';
end % END FOR

test = test./10e2; % Reduces the magnitude
test(:,2) = -test(:,1);

R = corrcoef(test);

%%
% EOF

laserEvents = ();

TrialRaw();
CorrCoef();
Waveform();
Rastergram();
PSTH()
SNR()
GenHeatmep();
FFT()
Spectrum()
Spectrogram()




%% Spectrum and Spectrogram

Fs = 1000;                    % Sampling frequency
T = 1/Fs;                     % Sample time
L = 1000;                     % Length of signal
t = (0:L-1)*T;                % Time vector
% Sum of a 50 Hz sinusoid and a 120 Hz sinusoid
x = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t); 
y = x + 2*randn(size(t));     % Sinusoids plus noise


NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

plot(f,2*abs(Y(1:NFFT/2+1))) 
