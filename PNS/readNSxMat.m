function DataStruct = readNSxMat(varargin)

% This function outputs a matrix of data by channel (row) from an nsxmat
% file (created using nsx2mat). The function takes up to two inputs. One of
% the inputs is the full filename to the NSx file to be converted. The
% other is a vector of channel numbers that specifies which channels to
% load. If this input is not provided, all channels in the data file will
% be loaded. The data is converted to class 'double', so be careful when
% loading all channels in the data file.
% 
% Example: DataStruct = readNSxMat('I:\sample.ns5',[1:4,9,64]);
% Example: DataStruct = readNSxMat([1,5,6]);
% Example: DataStruct = readNSxMat('I:\sample.ns5');
% Example: DataStruct = readNSxMat;
%
% Version Date: 20120717
% Author: Tyler Davis

FileName = '';
Channels = []; %vector of channels to output

switch nargin    
    case 2
        if isnumeric(varargin{1})
            FileName = varargin{2};
            Channels = varargin{1};
        else
            FileName = varargin{1};
            Channels = varargin{2};
        end            
    case 1
        if isnumeric(varargin{1})
            Channels = varargin{1};
        else
            FileName = varargin{1};
        end       
end

if isempty(FileName)
    [filename,path] = uigetfile('I:\Data\*.ns*','Choose nsxmat file...');
    FileName = fullfile(path,filename);
end

% Constructing output
load(FileName,'-mat','Header')
DataStruct = struct('header',[],'data',[],'dataID',[]);
DataStruct.header = Header;

ChannelID = mat2cell(Header.ChannelID,ones(1,size(Header.ChannelID,1)));
ChannelLabel = cellfun(@deblank,mat2cell(Header.ChannelLabel,ones(1,size(Header.ChannelLabel,1)),size(Header.ChannelLabel,2)),'uniformoutput',false);
UnitConversion = double(Header.MaxAnlgVal)./double(Header.MaxDigVal);
Units = cellfun(@deblank,mat2cell(Header.Units,ones(1,size(Header.Units,1)),size(Header.Units,2)),'uniformoutput',false);

if isempty(Channels)
    DataStruct.data = zeros(Header.ChannelCount,Header.ChannelLengthSamples);
    DataStruct.dataID = [ChannelID,ChannelLabel,Units];
    for k=1:Header.ChannelCount
        clc, disp(k)
        load(FileName,'-mat',['C',num2str(ChannelID{k})])
        eval(['DataStruct.data(',num2str(k),',:) = double(C',num2str(ChannelID{k}),')*UnitConversion(',num2str(k),');']);
        clear(['C',num2str(ChannelID{k})])
    end
else
    [~,cIdx,~] = intersect([ChannelID{:}],Channels);
    DataStruct.data = zeros(length(cIdx),Header.ChannelLengthSamples);
    DataStruct.dataID = [ChannelID(cIdx),ChannelLabel(cIdx),Units(cIdx)];
    for k=1:length(cIdx)
        clc, disp(ChannelID{cIdx(k)})
        load(FileName,'-mat',['C',num2str(ChannelID{cIdx(k)})])
        eval(['DataStruct.data(',num2str(k),',:) = double(C',num2str(ChannelID{cIdx(k)}),')*UnitConversion(',num2str(cIdx(k)),');']);
        clear(['C',num2str(ChannelID{cIdx(k)})])
    end
end

