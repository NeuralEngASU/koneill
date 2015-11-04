function mapFile = loadCMP(fileName)

% Loads an Blackrock provided mapfile data into a structure. This mapfile 
% contains information on the graphical placement of each channel, its 
% corresponding electrode and the bank and pin the electrode is connected
% to and it is used with many plotting functions.
%
% Usage:
%
%   electrodeMap = loadCMP(fileName);
%
%   WHERE:
%
%   filename:       Is the input filename of the CMP. If filename does not
%                   exist then the user will be prompted to select the
%                   file.
%
%   electrodeMap:   electrodeMap is the output structure of this file.
%
%   Kian Torab
%   kian.torab@utah.edu
%   Department of Bioengineering
%   University of Utah
%   Version 1.0.2 - March 8, 2010

if ~exist('fileName', 'var')
    [fileHandle pathHandle] = uigetfile('Y:\Code\Shared Codes\General Functions and LUTs\LUTs\Arrays\Cerebus Mapfiles\*.cmp');
    if ~fileHandle; disp('No file was selected.'); mapFile = []; return; end;
    mapfileDataCell = importdata([pathHandle fileHandle], ' ', 200);
else
    mapfileDataCell = importdata(fileName, ' ', 200);
end
mapfileDataCell(1:14) = [];
mapfileDataCellParsed = regexp(mapfileDataCell, '\t', 'split');


for i = 1:size(mapfileDataCellParsed, 1)
    mapfileParsed(i).Column   = str2num(mapfileDataCellParsed{i,:}{1});
    mapfileParsed(i).Row      = str2num(mapfileDataCellParsed{i,:}{2});
    mapfileParsed(i).Bank     = mapfileDataCellParsed{i,:}{3}-'@';
    mapfileParsed(i).Pin      = str2num(mapfileDataCellParsed{i,:}{4});
    ElectNum = str2num(mapfileDataCellParsed{i,:}{5}(5:end));
    if isempty(ElectNum)
        mapfileParsed(i).ElectNum = str2num(mapfileDataCellParsed{i,:}{5}(2:3));
    else
        mapfileParsed(i).ElectNum = ElectNum;
    end
    mapfileParsed(i).ChanNum  = (mapfileParsed(i).Bank - 1) * 32 + mapfileParsed(i).Pin;
    mapfileParsed(i).Label    = mapfileDataCellParsed{i,:}{5};
end

mapFile([mapfileParsed.ChanNum]) = mapfileParsed;

if ~nargout
    assignin('caller', 'mapfile', mapFile);
end