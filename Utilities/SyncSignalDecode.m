%% Sync Signal Decode

function timeMat = SyncSignalDecode(syncData, Fs)
% Function setup
maxSyncTSLength = floor(4.3*Fs); % Units samples. Timestamps in sync signal are 43 bits long and sent at a rate of 10 Hz.
minSyncTSPeriod = floor(10*Fs); % Units samples. Timestamps in sync signal can have a minimum period of 10 sec.
% BitLength = floor(0.1*Fs); % Units samples.
bitLength = floor(maxSyncTSLength/43);

% Finding voltages (Values) for each timestamp (TS) in the SyncSignal
params.Fs = Fs;
params.eventGap = 3; % [seconds] Gap threshold between event groups
eventIdx = FindStepEvents(syncData, params);

% Compute timestamp
strMat = [];
% loop over each event
for ii = 1:size(eventIdx,2)-1
    testEvent = eventIdx{ii};
    testEvent = testEvent-testEvent(1)+1;
    if size(testEvent,2) <= 2
        timeStamp{ii} = nan;
    else
        signal = zeros(1, maxSyncTSLength);
        for jj = 1:2:size(testEvent,2)
            if jj+1>size(testEvent,2)
                signal(testEvent(jj):end) = 1;
            else
                signal(testEvent(jj):testEvent(jj+1)) = 1;
            end
        end
        
        signal2 = reshape(signal(1:bitLength*43), bitLength, 43);
        syncBits = round(mean(signal2,1));
%         if syncBits(end) == 0
%             disp(ii)
%             syncBits = ~syncBits;
%         end
        syncBitsStr = num2str(fliplr(syncBits(:,12:end))')';
        syncBitsStr = syncBitsStr(1:3:size(syncBitsStr,1),:);
        syncSec = bin2dec(syncBitsStr);
        
        strMat = [strMat;syncBitsStr];
        
        % Convert from Labview to Matlab time (64800 represents seconds in 18 hours)
        labview2MatlabOffset = etime(datevec('01-01-1904 00:00:00','mm-dd-yyyy HH:MM:SS'),datevec('01-01-0000 00:00:00','mm-dd-yyyy HH:MM:SS')) + 64800;

        labViewTime = datenum('19040101 00:00:00', 'yyyymmdd HH:MM:SS');
        expTime = (labview2MatlabOffset + syncSec)/8.64e4;
        matlabTime = datenum('0000101 00:00:00', 'yyyymmdd HH:MM:SS');
        offset = labViewTime - matlabTime;
        
%         SyncTSDateVec(ii,:) = datevec(expTime);
        
        timeMat(ii) = expTime;

        timeStamp{ii} = expTime;
    end   
end % END FOR

% for kk = 1:size(timeStamp,2)
%     if ~isnan(timeStamp{kk})
%         disp(datestr(timeStamp{kk}, 'yyyymmdd HH:MM:SS'))
%     end
% end % END FOR

end % END FUNCTION SyncSignalDecode

% EOF