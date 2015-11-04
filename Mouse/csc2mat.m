function [ output_args ] = csc2mat( input_args )

NSxMatFullName = [fullfile(NSxPath,NSxName),NSxExt,'mat'];
TMPFullName = [fullfile(NSxPath,NSxName),'.tmp'];

BegOfData = ftell(FID);
fseek(FID, 0, 'eof');
EndOfData = ftell(FID);
fseek(FID, BegOfData, 'bof');

Header.DataLengthBytes = EndOfData - BegOfData;
Header.ChannelLengthBytes = Header.DataLengthBytes/Header.ChannelCount;
Header.ChannelLengthSamples = Header.DataLengthBytes/2/Header.ChannelCount;

% Determining system memory to maximize data segments
SystemMemory = regexp(evalc('feature memstats'),'\d*(?= MB)','match');
SystemMemory = str2double(SystemMemory{2})*1e6; % Units bytes

% Calculating maximum data segment to load into memory
SegmentCount = ceil(Header.DataLengthBytes/(0.75*SystemMemory));        
SegmentSamples = round(Header.ChannelLengthSamples/SegmentCount);
SegmentDivisor = floor(Header.ChannelLengthSamples/SegmentSamples);
SegmentRemainder = rem(Header.ChannelLengthSamples,SegmentSamples);
if SegmentRemainder==0
    SegmentMatrix = repmat(SegmentSamples,SegmentDivisor,1);
else
    SegmentMatrix = [repmat(SegmentSamples,SegmentDivisor,1);SegmentRemainder];
end

% Parsing and saving to *.tmp file
save(TMPFullName,'Header','-v7.3')
for k = 1:size(SegmentMatrix,1) 
    clc, fprintf('NSX2MAT Parsing: %0.0f%% complete\nLoading %0.1f GB of data\n',((k-1)*Header.ChannelCount*100)/(size(SegmentMatrix,1)*Header.ChannelCount),SegmentMatrix(k)*2/1e7)
    tempData = fread(FID, [Header.ChannelCount,SegmentMatrix(k)], '*int16');    
    for m = 1:Header.ChannelCount 
        clc, fprintf('NSX2MAT Parsing: %0.0f%% complete\n',(((k-1)*Header.ChannelCount+m)*100)/(size(SegmentMatrix,1)*Header.ChannelCount)) 
        tempSubChan = ['C',num2str(Header.ChannelID(m)),'_',num2str(k)];
        eval([tempSubChan,' = tempData(m,:);']);
        save(TMPFullName,tempSubChan,'-append','-v7.3')                                   
        clear(tempSubChan)
    end
    clear('tempData')
end



end

