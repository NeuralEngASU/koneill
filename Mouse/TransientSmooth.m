function [ smoothData ] = TransientSmooth( inputData )

% Get size of data
[r,c] = size(inputData);

% Initialize NaN matrix
smoothData = nan(r, c);

% Manually compute the first three points
smoothData(1,:) = inputData(1,:);
smoothData(2,:) = nanmean(inputData(1:3, :));
smoothData(3,:) = nanmean(inputData(1:5, :));

% Average over n-2, n-1, n, n+1, n+2 and save data to n
for i = 4:r-2

    smoothData(i,:) = nanmean(inputData((i-2):(i+2), :));

end % END FOR

% manually compute the last 2 points
smoothData(end-1,:) = nanmean(inputData(end-2:end, :));
smoothData(end,:) = inputData(end,:);

end % END FUNCTION

% EOF