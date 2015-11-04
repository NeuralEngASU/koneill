function nsx2mat(varargin)

% Converts NSx files to MAT files and saves as *.nsxmat within the
% specified directory. The function takes up to two inputs. One of the
% inputs is the full filename to the NSx file to be converted. Local disk
% space triple the size of the file is required for this function to run.
% Works for filespec 2.2 and 2.3. If you want to downsample, add the
% downsample frequency (i.e. '2000') in units S/sec as a string as an
% additional input. 
%
% Example: nsx2mat('I:\sample.ns5','10000');
% Example: nsx2mat('2000','I:\sample.ns4');
% Example: nsx2mat('I:\sample.ns4');
% Example: nsx2mat('1000');
% Example: nsx2mat;
%
% Version Date: 20120808
% Version comments: Added the ability to read split files and removed filespec 2.1
% Author: Tyler Davis

FileName = '';
dsFs = []; %downsample frequency
nsx2mat_version = 20120808;

switch nargin    
    case 2
        if isnumeric(varargin{1})
            FileName = varargin{2};
            dsFs = varargin{1};
        else
            FileName = varargin{1};
            dsFs = varargin{2};
        end            
    case 1
        if isnumeric(varargin{1})
            dsFs = varargin{1};
        else
            FileName = varargin{1};
        end       
end

if isempty(FileName)
    [filename,path] = uigetfile('I:\Data\*.ns*','Choose nsx file...');
    FileName = fullfile(path,filename);
end

NSxFullName = FileName;
NSxMatFullName = [NSxFullName,'mat'];
TMPFullName = [NSxFullName(1:end-3),'tmp'];

switch NSxFullName(end-2:end)
    case 'ns5'
        Fs = 30000;
    case 'ns4'
        Fs = 10000;
    case 'ns3'
        Fs = 2000;
    case 'ns2'
        Fs = 1000;
    otherwise
        disp('Choose an NSx file')
        return
end

% Version and sample frequency checking
try
    if exist(NSxMatFullName,'file')
        load(NSxMatFullName,'-mat','Header')
        if isfield(Header,'nsx2mat_version')
            if Header.nsx2mat_version~=nsx2mat_version
                delete(NSxMatFullName);
            elseif ~isempty(dsFs)
                if Header.Fs~=dsFs
                    delete(NSxMatFullName);
                end
            else
                if Header.Fs~=Fs
                    delete(NSxMatFullName);
                end
            end
        else
            delete(NSxMatFullName);
        end
        clear('Header')
    end
catch ME
    disp('Cannot delete old nsxmat file. Run Matlab as administrator.')
    return
end

% Checking if nsxmat file already exists
if exist(NSxMatFullName,'file')
    disp('nsxmat file already exists!')
    return
end

% Checking filespec
FID = fopen(NSxFullName, 'r', 'l');
fseek(FID, 8, 'bof');
filespec = fread(FID, [1,2], '*uchar');
fseek(FID, 0, 'bof');
if ~(all(filespec==[2,2]) || all(filespec==[2,3]))
    disp('Can not read this filespec!')
    return
end

% Creating filter object
if ~isempty(dsFs)
    dsFs = Fs/round(Fs/dsFs);
    if Fs > dsFs
        N = 10; %order
        Fpass = dsFs/4; %passband frequency
        Apass = 1; %passband ripple (dB)
        Astop = 80; %stopband attenuation (dB)
        h = fdesign.lowpass('N,Fp,Ap,Ast',N,Fpass,Apass,Astop,Fs);
        Hd = design(h,'ellip');
    else
        disp('Cannot downsample to this frequency. Choose a different value.')
        return
    end
end

Header.nsx2mat_version = nsx2mat_version;

% Reading NSx file header
Header.FileID       = fread(FID, [1,8],   '*char');
Header.FileSpec     = fread(FID, [1,2],   '*uchar');
Header.HeaderBytes  = fread(FID, [1,1],   '*uint32');
Header.Fs           = fread(FID, [1,16],  '*char');
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

BegOfDataHeader = ftell(FID);
fseek(FID, 0, 'eof');
EndOfFile = ftell(FID);
fseek(FID, BegOfDataHeader, 'bof');

% Checking for multiple data headers
k = 1;
while ftell(FID)~=EndOfFile
    DataHeader(k,1) = fread(FID, [1,1], '*uint8');
    if DataHeader(k,1)~=1
        disp('Error reading data headers!')
        return
    end
    DataTimestamp(k,1) = fread(FID, [1,1], '*uint32');
    ChannelLengthSamples(k,1) = fread(FID, [1,1], 'uint32=>double');
    BegOfData(k,1) = ftell(FID);
    fseek(FID,ChannelLengthSamples(k,1)*Header.ChannelCount*2,'cof');
    DataLengthBytes(k,1) = ftell(FID) - BegOfData(k,1);
    k = k+1;
end

% Check if pauses exist in data
if (length(DataHeader)==2 && ChannelLengthSamples(1)~=1) || length(DataHeader)>2
    disp('Pauses exist in this data set! This version of nsx2mat cannot parse paused data')
    return
end 

% Check if data length in header is equal to calculated data length
if DataLengthBytes(end)~=ChannelLengthSamples(end)*Header.ChannelCount*2
    disp('Header and calculated data lengths are different!')
    return
end

% Move back to beginning of last data segment
fseek(FID,BegOfData(end),'bof'); 

% Discarding info from extra data header and updating main header
% An extra data header and a single sample of data are found when files are automatically split using firmware 6.03
Header.ChannelLengthSamples = ChannelLengthSamples(end);

% Determining system memory to maximize data segments
SystemMemory = regexp(evalc('feature memstats'),'\d*(?= MB)','match');
SystemMemory = str2double(SystemMemory{2})*1e6; % Units bytes

% Calculating maximum channel samples to load into memory
SegmentCount = ceil(DataLengthBytes(end)/(0.75*SystemMemory));
SegmentMatrix = repmat(floor(Header.ChannelLengthSamples/SegmentCount),SegmentCount,1);
SegmentMatrix(1) = SegmentMatrix(1) + rem(Header.ChannelLengthSamples,SegmentCount);

% Parsing and saving to *.tmp file
save(TMPFullName,'Header','-v7.3')
for k = 1:size(SegmentMatrix,1)
    clc, fprintf('NSX2MAT Parsing: %0.0f%% complete\nLoading %0.1f GB of data\n',((k-1)*Header.ChannelCount*100)/(size(SegmentMatrix,1)*Header.ChannelCount),SegmentMatrix(k)*2*Header.ChannelCount/1e9)
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

% Modifying header for split files
if DataTimestamp(end)>0
    Header.ChannelLengthSamples = Header.ChannelLengthSamples + DataTimestamp(end) - 1;
end

% Modifying header for downsampled files
if ~isempty(dsFs)
    Header.ChannelLengthSamples = (Header.ChannelLengthSamples-rem(double(Header.ChannelLengthSamples)-1,round(Fs/dsFs))-1)/round(Fs/dsFs)+1;
    Header.Fs = dsFs;
    Header.Period = round(double(Header.Resolution)/dsFs);
else
    Header.Fs = Fs;
end

% Reconstructing from *.tmp file
save(NSxMatFullName,'Header','-v7.3')
for m = 1:Header.ChannelCount
    clc, fprintf('NSX2MAT Reconstructing: %0.0f%% complete\n',m*100/double(Header.ChannelCount))
    tempChan = ['C',num2str(Header.ChannelID(m))];
    eval([tempChan,' = [];'])
    for k = 1:size(SegmentMatrix,1)
        tempSubChan = ['C',num2str(Header.ChannelID(m)),'_',num2str(k)];
        load(TMPFullName,tempSubChan,'-mat')
        eval([tempChan,' = [',tempChan,',',tempSubChan,'];'])
        clear(tempSubChan)
    end
    %%%%%%%%%%%% Buffer with zeros for split files %%%%%%%%%%%%%%%%%%%%
    if DataTimestamp(end)>0
        eval([tempChan,' = [zeros(1,DataTimestamp(end)-1),',tempChan,'];'])
    end
    %%%%%%%%%%%%%%%%%%%%%%%% Downsample %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(dsFs)
        eval([tempChan,' = filtfilt(Hd.sosMatrix,Hd.ScaleValues,double(',tempChan,'));'])
        eval([tempChan,' = ',tempChan,'(1:',num2str(round(Fs/dsFs)),':end);'])
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save(NSxMatFullName,tempChan,'-append','-v7.3')
    clear(tempChan)
end

delete(TMPFullName)
fclose(FID);


