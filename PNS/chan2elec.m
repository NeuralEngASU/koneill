function Electrode = chan2elec(Channel,cmpFile)

% Version Date: 9-22-2010
% Author: Tyler Davis

if nargin==1
    mapFile = loadCMP('DefaultArrayMap.cmp');
else
    mapFile = loadCMP(cmpFile);
end

ChanElecLUT = [[mapFile.ChanNum]',[mapFile.ElectNum]'];
ChanElecLUT(97:100,1) = 97:100;
ChanElecLUT(97:100,2) = setdiff(1:100,ChanElecLUT(:,2));

Electrode = ChanElecLUT(Channel,2);