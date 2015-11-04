function [ slopeData ] = TransientSlope( inputData )

% Get the size of the date
[r,c] = size(inputData);

% Create NaN matrix
slopeData = nan(r,c);

for i = 1:r-1

    slopeData(i,:) = inputData(i+1, :) - inputData(i,:);
    
end % END FOR

end % END FUNCTION

% EOF