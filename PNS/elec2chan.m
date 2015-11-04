function Channel = elec2chan(Electrode,cmpFile)

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
ChanElecLUTsorted = sortrows(ChanElecLUT,2);

Channel = ChanElecLUTsorted(Electrode,1);