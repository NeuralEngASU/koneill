%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLI_Test
%   Author: Kevin O'Neill
%   Date: 2015.03.16
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load Data %%

% Make sure the data files are in MATLAB's path.

NEV = openNEV('E:\data\PLI\201203\20120627-164445\20120627-164445-001.nev');

nsx2mat('E:\data\PLI\201203\20120627-164445\20120627-164445-001.ns4')
load('E:\data\PLI\201203\20120627-164445\20120627-164445-001.ns4mat', '-mat')
%% Phase Lag Index
% This calculates the PLI between two signal sources. The calculation is
% the difference in phase between two sources. If the PLI is 0, there is no
% coupling, or coupling around 0 mod pi. If PLI is one there is perfect
% phase locking at a value of delta(phi) from 0 mod pi. 0 < PLI < 1

tIdx = 1:100000;
sig1 = sin(pi*tIdx/100);
sig2 = sin(pi*tIdx/100);

sig1FFT = fft(sig1);
sig2FFT = fft(sig2);

sig1FFT = sig1FFT(2:end);
sig2FFT = sig2FFT(2:end);

% sig1FFT(logical(abs(sig1FFT < 1e-10))) = 0;
% sig2FFT(logical(abs(sig2FFT < 1e-10))) = 0;

% deltaPhi = angle(sig1FFT) - angle(sig2FFT);

% sig1FFT = fft(double(C1(tIdx)));
% sig2FFT = fft(double(C25(tIdx)));
% 
% sig1FFT = sig1FFT(2:end);
% sig2FFT = sig2FFT(2:end);

% sig1FFT(logical(abs(sig1FFT))) = 0;
% sig2FFT(logical(abs(sig2FFT))) = 0;
deltaPhi = angle(sig1FFT) - angle(sig2FFT);

PLI = abs(mean(sign(deltaPhi)));

% plot(C1(tIdx)); hold on; plot(C2(tIdx)); hold off

%% Hilbert method


tIdx = 1:1000;%000;
% sig1 = sin(pi*tIdx/100);
% sig2 = sin(pi*tIdx/100);


% phaseAngle =  linspace(0,pi,10000);

% for kk = 1:10000
kk = 1;
% sig1 = sin(pi*tIdx/100);
% sig2 = sin(pi*tIdx/100 + phaseAngle(kk));

sig1 = double(C1(tIdx));
sig2 = double(C1(tIdx));

% sig1Phi = angle(atan(hilbert(sig1)./sig1));
% sig2Phi = angle(atan(hilbert(sig2)./sig2));

sig1Phi = atan(imag(hilbert(sig1))./sig1);
sig2Phi = atan(imag(hilbert(sig2))./sig2);

deltaPhi = sig1Phi - sig2Phi;

testR(kk) = abs(1/length(deltaPhi) * sum( exp(1i * deltaPhi(1:end-1))) );

PLI(kk) = abs(mean(sign(deltaPhi)));

% end % END FOR
% figure(3)
% plot(sig1./max(sig1)); hold on; plot(sig2./max(sig2)); %plot(deltaPhi./max(deltaPhi));hold off

% deltaPhi(deltaPhi > pi) = deltaPhi(deltaPhi > pi) - 2*pi; 

%% Surrogate Data

tIdx = 1:10000;
sig1 = double(C1(tIdx));
sig2 = double(C4(tIdx));

surrDataCount = 100;

surrData = zeros(surrDataCount, length(sig1));
phaseDiff = -pi + (pi+pi)*rand(1,surrDataCount);

R = zeros(1,surrDataCount);
PLI = zeros(1,surrDataCount);

sig1Phi = atan(imag(hilbert(sig1))./sig1);
fs = length(tIdx);

for jj = 1:surrDataCount
    
    surrData(jj,:) = fft(sig2);
    f = fftshift(( 0:(fs-1)) - fs/2);
    surrData(jj,:) = surrData(jj,:).* exp(1i*2*f*phaseDiff(jj));   
    surrData(jj,:) = real(ifft(surrData(jj,:)));
    
    
    sig2Phi = atan(imag(hilbert(surrData(jj,:)))./surrData(jj,:));

    deltaPhi = sig1Phi - sig2Phi;
    
    R(jj) = abs(1/length(deltaPhi) * sum( exp(1i * deltaPhi(1:end-1))) );
    PLI(jj) = abs(mean(sign(deltaPhi)));
    
end % END IF

sig2Phi = atan(imag(hilbert(sig2))./sig2);
    
deltaPhi = sig1Phi - sig2Phi;

R(end+1) = abs(1/length(deltaPhi) * sum( exp(1i * deltaPhi(1:end-1))) );
PLI(end+1) = abs(mean(sign(deltaPhi)));

zscore((PLI - PLI(end))')

%% Surrogate Data Full Length

% Allocate variables
sig1 = zeros(1,length(C1));
sig2 = zeros(1,length(C1));

numChannels = 32;
spanChannels = 1:numChannels;

surrDataCount = 20;

surrData = zeros(1, length(sig1));
phaseDiff = -pi + (pi+pi)*rand(numChannels,numChannels,surrDataCount);

R = zeros(numChannels,numChannels,surrDataCount+1);
PLI = zeros(numChannels,numChannels,surrDataCount+1);

fs = length(sig1);

zScoreReport = zeros(numChannels,numChannels,surrDataCount+1);

for ii = spanChannels
    for jj = spanChannels
        
        eval(['sig1 = double(C', num2str(ii),');']);
        eval(['sig2 = double(C', num2str(jj),');']);
        
        sig1Phi = atan(imag(hilbert(sig1))./sig1);
        
        for kk = 1:surrDataCount

            tic;
            
            surrData = fft(sig2);
            f = fftshift(( 0:(fs-1)) - fs/2);
            surrData = surrData.* exp(1i*2*f*phaseDiff(ii,jj,kk));
            surrData = real(ifft(surrData));
                        
            sig2Phi = atan(imag(hilbert(surrData))./surrData);
            
            deltaPhi = sig1Phi - sig2Phi;
            
            R(ii,jj,kk) = abs(1/length(deltaPhi) * sum( exp(1i * deltaPhi(1:end-1))) );
            PLI(ii,jj,kk) = abs(mean(sign(deltaPhi)));
            timeSpent = toc;
            
            numTotal =  surrDataCount*length(spanChannels)^2;
            numComplete = (ii-1)*numChannels*surrDataCount + (jj-1)*surrDataCount + (kk);
            timeLeft = (numTotal - numComplete) * timeSpent;
            
            timeSpent = numComplete * timeSpent;
            
            clc;
            disp(sprintf('%d/%d', numComplete, numTotal));
            disp(sprintf('Time left: %f, minutes',timeLeft/60))
            disp(sprintf('Time spent: %f, minutes',timeSpent/60))
            
        end % END FOR surrDataCount  
        
        sig2Phi = atan(imag(hilbert(sig2))./sig2);
        
        deltaPhi = sig1Phi - sig2Phi;
               
        R(ii,jj,end) = abs(1/length(deltaPhi) * sum( exp(1i * deltaPhi(1:end-1))) );
        PLI(ii,jj,end) = abs(mean(sign(deltaPhi)));
        
        zScoreReport(ii,jj,:) = zscore(PLI(ii,jj,:) - PLI(ii,jj,end));
        
    end % END FOR sig2 spanChannels   
end % END FOR sig1 spanChannels

% Flip the third dimension in order to bring the real value to the first
% layer
% R            = permute(           R, [1,2,-3]);
% PLI          = permute(         PLI, [1,2,-3]);
% zScoreReport = permute(zScoreReport, [1,2,-3]);

R2            = R(:,:,end:-1:1);
PLI2          = PLI(:,:,end:-1:1);
zScoreReport2 = zScoreReport(:,:,end:-1:1);


disp('Done')
% EOF