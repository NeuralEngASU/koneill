% function [ output_args ] = PLISeparateBands( data, class )

numSamp = size(data,1);
numChan = size(data,2);
numTrial = size(data,3);

bands = [ 1,   4;... % Delta Band  
          4,   8;... % Theta
          8,  12;... % Alpha
         12,  30;... % Beta
         30,  80;... % Gamma
         80, 200];   % Chi

bandNames = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma', 'Chi'};

parpool(16)

Apass =  1; % dB, Passband attenuation
Astop1 = 60; % dB, Stopband attenuation
Astop2 = 60; % dB, Stopband attenuation
order = 20; % The order of the filter
match  = 'both';  % Band to match exactly

% bandData = zeros(size(data,1), size(data,2), size(data,3));

Header.bands = bands;
Header.bandNames = bandNames;
Header.class = class;
Header.Fs = 5000;
Header.filtType = 'butter';
Header.stopBand = 0.9;
%%
timeClock = tic;

for b = 1:size(bands,1)
    bandData = zeros(size(data,1), size(data,2), size(data,3));

    Header.currentBand = b;
    
    Fstop1 = bands(b,1)-0.9;
    Fstop2 = bands(b,2)+0.9;
    
    Fpass1 = bands(b,1);
    Fpass2 = bands(b,2);
    
    for t = 1:numTrial
        
        % Construct an FDESIGN object and call its ELLIP method.
        h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, Astop2, Header.Fs);
        Hd = design(h, 'ellip', 'MatchExactly', match);
%         d = fdesign.bandpass('N,F3dB1,F3dB2',60,Fstop1,Fstop2,5000);
%         Hd = design(d,'ellip'); %fvtool(Hd)
        parfor c = 1:numChan
            h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, Astop2, Header.Fs);
            Hd = design(h, 'ellip', 'MatchExactly', match);
%             d = fdesign.bandpass('N,F3dB1,F3dB2',60,Fstop1,Fstop2,5000);
%             Hd = design(d,'ellip'); %fvtool(Hd)
            % h=fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, Astop2);
            % Hd=design(h, 'butter');
            bandData(:,c,t) = filtfilt(Hd.sosMatrix,Hd.ScaleValues,data(:,c,t));
        end
        timeElap = toc(timeClock);
        
        timeTrial = timeElap/(t + numTrial*(b-1));
        
        clc
        fprintf('Current Band: %d/%d, %s\n', b, size(bands,1), bandNames{b});
        fprintf('Current Trial: %d/%d\n', t, numTrial);
        fprintf('Time Per Trial: %0.3f seconds\n',timeTrial);
        fprintf('Time Spent: %0.3f minutes\n', timeElap / 60);
        fprintf('Time Left: %0.3f minutes\n', timeTrial * ((6-b)*numTrial + numTrial-t) / 60 );
        
    end % END FOR
    
    bandData([1:round(0.5*Header.Fs), end-round(0.5*Header.Fs)+1:end],:,:) = [];
    
    save(['E:\data\PLI\delta\PLIOutput\SeperateBands\Delta_PLI_', bandNames{b}, '.mat'], 'bandData', 'Header', '-v7.3');
end %END FOR
     
delete(gcp)
% end % END FUNCTION

% EOF