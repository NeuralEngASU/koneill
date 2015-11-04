function [ patchCoord ] = GenPatchCoord1( numChans )

patchCoord = zeros(4, numChans, 2);

xVals = zeros(4, numChans, 1);
yVals = zeros(4, numChans, 1);

plotHeight = [numChans:-1:0];
plotWidth  = [0:1:numChans];
for i = 1:numChans

    xVals(1,i) = plotWidth(i);
    xVals(2,i) = plotWidth(i+1);
    xVals(3,i) = plotWidth(i+1);
    xVals(4,i) = plotWidth(i);
    
    yVals(1,i) = plotHeight(i);
    yVals(2,i) = plotHeight(i);
    yVals(3,i) = plotHeight(i+1);
    yVals(4,i) = plotHeight(i+1);
    
end % END FOR

patchCoord(:,:,1) = xVals;
patchCoord(:,:,2) = yVals;

end % END FUNCTION