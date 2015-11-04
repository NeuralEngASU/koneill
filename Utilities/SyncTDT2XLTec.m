function timeStampTDT = SyncTDT2XLTec(syncDataTDT, Fs)

% Extract only rises
syncData = zeros(1, size(syncDataTDT,2));

syncData(diff(syncDataTDT)==1) = 1;

% Generate Times
timeStampTDT = SyncSignalDecode(syncData, Fs);

end % END FUNCTION SyncTDT2XLTec

% EOF