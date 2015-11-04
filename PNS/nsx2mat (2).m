function nsx2mat(varargin)

% Converts NSx files to MAT files and saves as *.nsxmat within the
% specified directory. The function takes up to two inputs. One of the
% inputs is the full filename to the NSx file to be converted. Local disk
% space triple the size of the file is required for this function to run.
% Works for filespec 2.1 and 2.2. If you want to downsample an ns5 to
% 2kS/sec, add the string '-ds' as a second input. If not, leave this input
% blank.
%
% Example: nsx2mat('I:\sample.ns5','-ds');
% Example: nsx2mat('I:\sample.nsx');
% Example: nsx2mat('-ds');
% Example: nsx2mat('I:\sample.nsx');
% Example: nsx2mat;
%
% Version Date: 20111212
% Author: Tyler Davis

NSxFullName = '';
downsampleFlag = '';

switch nargin
    case 2
        if exist(varargin{1},'file')==2
            NSxFullName = varargin{1};
            if strcmp(varargin{2},'-ds')
                downsampleFlag = varargin{2};            
            end
        elseif strcmp(varargin{1},'-ds')
            downsampleFlag = varargin{1};
            if exist(varargin{2},'file')==2
                NSxFullName = varargin{2};            
            end
        else
            if exist(varargin{2},'file')==2
                NSxFullName = varargin{2};
            elseif strcmp(varargin{2},'-ds')
                downsampleFlag = varargin{2};            
            end
        end        
    case 1
        if exist(varargin{1},'file')==2
            NSxFullName = varargin{1};
        elseif strcmp(varargin{1},'-ds')
            downsampleFlag = varargin{1};        
        end    
end

if strcmp(downsampleFlag,'-ds')
    % Elliptic Lowpass filter. All frequency values are in Hz.
    Fs    = 30000;% Sampling Frequency    
    N     = 10;   % Order
    Fpass = 500;  % Passband Frequency
    Apass = 1;    % Passband Ripple (dB)
    Astop = 80;   % Stopband Attenuation (dB)    
    h = fdesign.lowpass('N,Fp,Ap,Ast', N, Fpass, Apass, Astop, Fs);
    Hd = design(h, 'ellip');   
end

if isempty(NSxFullName)
    [filename,path] = uigetfile('I:\Data\*.ns*','Choose NSx file...');
    NSxFullName = fullfile(path,filename);
end

[NSxPath,NSxName,NSxExt] = fileparts(NSxFullName);
NSxMatFullName = [fullfile(NSxPath,NSxName),NSxExt,'mat'];
TMPFullName = [fullfile(NSxPath,NSxName),'.tmp'];

% Checking filespec
FID = fopen(NSxFullName, 'r', 'l');
fseek(FID, 8, 'bof');
filespec = fread(FID, [1,2],   '*uchar');
fseek(FID, 0, 'bof');

% Reading NSx file
if exist(NSxMatFullName,'file')==0
    
    if all(filespec==[2,2])
        
        % Reading filespec 2.2
        Header.FileID       = fread(FID, [1,8],   '*char');
        Header.FileSpec     = fread(FID, [1,2],   '*uchar');
        Header.HeaderBytes  = fread(FID, [1,1],   '*uint32');
        Header.Label        = fread(FID, [1,16],  '*char');
        Header.Comment      = fread(FID, [1,256], '*char');
        Header.Period       = fread(FID, [1,1],   '*uint32');
        Header.Resolution   = fread(FID, [1,1],   '*uint32');
        Header.TimeOrigin   = fread(FID, [1,8],   '*uint16');
        Header.ChannelCount = fread(FID, [1,1],   'uint32=>double');
        
        for k = 1:Header.ChannelCount
            Header.Type(k,:)           = fread(FID, [1,2],  '*char');
            Header.ChannelID(k,:)      = fread(FID, [1,1],  '*uint16');
            Header.ChannelLabel(k,:)   = fread(FID, [1,16], '*char');
            Header.PhysConnector(k,:)  = fread(FID, [1,1],  '*uint8');
            Header.ConnectorPin(k,:)   = fread(FID, [1,1],  '*uint8');
            Header.MinDigVal(k,:)      = fread(FID, [1,1],  '*int16');
            Header.MaxDigVal(k,:)      = fread(FID, [1,1],  '*int16');
            Header.MinAnlgVal(k,:)     = fread(FID, [1,1],  '*int16');
            Header.MaxAnlgVal(k,:)     = fread(FID, [1,1],  '*int16');
            Header.Units(k,:)          = fread(FID, [1,16], '*char');
            Header.HighFreqCorner(k,:) = fread(FID, [1,1],  '*uint32');
            Header.HighFreqOrder(k,:)  = fread(FID, [1,1],  '*uint32');
            Header.HighFiltType(k,:)   = fread(FID, [1,1],  '*uint16');
            Header.LowFreqCorner(k,:)  = fread(FID, [1,1],  '*uint32');
            Header.LowFreqOrder(k,:)   = fread(FID, [1,1],  '*uint32');
            Header.LowFiltType(k,:)    = fread(FID, [1,1],  '*uint16');
        end
        
        Header.DataHeader     = fread(FID, [1,1], '*uint8');
        Header.DataTimestamp  = fread(FID, [1,1], '*uint32');
        Header.DataLength     = fread(FID, [1,1], '*uint32');
        
    else
        
        % Reading filespec 2.1
        Header.FileID       = fread(FID, [1,8],  '*char');
        Header.Label        = fread(FID, [1,16], '*char');
        Header.Period       = fread(FID, [1,1],  '*uint32');
        Header.ChannelCount = fread(FID, [1,1],  'uint32=>double');
        Header.ChannelID    = fread(FID, [Header.ChannelCount,1], '*uint32');
        
    end

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

% Reconstructing from *.tmp file
if strcmp(downsampleFlag,'-ds')
    Header.Downsampled = '2kS/s';   
end
save(NSxMatFullName,'Header','-v7.3')
for m = 1:Header.ChannelCount
    clc, fprintf('NSX2MAT Reconstructing: %0.0f%% complete\n',m*100/Header.ChannelCount)
    tempChan = ['C',num2str(Header.ChannelID(m))];
    eval([tempChan,' = [];'])
    for k = 1:size(SegmentMatrix,1)
        tempSubChan = ['C',num2str(Header.ChannelID(m)),'_',num2str(k)];
        load(TMPFullName,tempSubChan,'-mat')
        eval([tempChan,' = [',tempChan,',',tempSubChan,'];'])
        clear(tempSubChan)
    end
    %%%%%%%%%%%%%%%%%%%%%%%% Downsample %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(downsampleFlag,'-ds')
        % Filtering data with 500Hz lowpass elliptical and downsampling to 2kS/sec
        eval([tempChan,' = filtfilt(Hd.sosMatrix,Hd.ScaleValues,double(',tempChan,'));'])
        eval([tempChan,' = ',tempChan,'(1:15:end);'])
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save(NSxMatFullName,tempChan,'-append','-v7.3')
    clear(tempChan)
end

delete(TMPFullName)
fclose(FID);

else
    clc, fprintf('nsxmat file already exists!\n')
end


