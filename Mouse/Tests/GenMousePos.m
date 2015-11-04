function [ varargout ] = GenMousePos( mouseData, varargin )

mouseData = round(mouseData); % Rounds the mouse data to the nearest integer

[row, col] = size(mouseData); % Finds the size of mouseData. Should be Nx2

if col ~= 2
    fprintf(['The given data is not of the right size.\n', ...
        'The matrix should be Nx2 where column 1 contains the x-values and column 2 contains the y-values.\n']);
    return;
end % END IF

% 1/2 the size of the gaussian. +/- sizeGauss
sizeGauss = varargin{1};

x = linspace(-sizeGauss, sizeGauss, 2*sizeGauss+1);
y = x';               
[X,Y] = meshgrid(x, y);
z = exp(-( X.^2 + Y.^2) / (2 * 1000));
% surf(x,y,z);shading interp

% Initializes the matrix to contain the heatmap data.
% This matrix needs a buffer zone so the size of the gaussian
% won't cause an error.
heatValMat = zeros(1080 + 2*sizeGauss+1, 1920 + 2*sizeGauss + 1);  

% Loops over every row of mouseData
for i = 1:size(mouseData, 1)
    
    % For every iteration, add the size of the gaussian buffer+1
    xPos = mouseData(i,1) + sizeGauss+1;
    yPos = mouseData(i,2) + sizeGauss+1;
    
    % Find the bouding box around the x-y-position.
    yMat = [yPos-sizeGauss:yPos+sizeGauss];
    xMat = [xPos-sizeGauss:xPos+sizeGauss];
    
    % Add the z function values to the bounding box
    heatValMat(yMat, xMat) = heatValMat(yMat, xMat) + z;
    
end % END FOR

% Calculate the total area added to the heatValMat per sample (valStep)
% and calculate the total area added to heatValMat (valTotal)
valStep = sum(z(:));
valTotal = sum(heatValMat(:));
% valTotal/valStep should equal the number of rows in mouseData

% Calculates the total "time" the mouse was in each section of the cage
lengthMat = size(heatValMat, 2);
Fs = 20;

timeSect(1) = sum(sum(heatValMat(:, 1                        : floor(lengthMat/3)   )));
timeSect(2) = sum(sum(heatValMat(:, floor(lengthMat/3) + 1   : floor(lengthMat/3)*2 )));
timeSect(3) = sum(sum(heatValMat(:, floor(lengthMat/3)*2 + 1 : end                  )));

timeSect = timeSect./valStep;
timeSect = timeSect./Fs;

timeSect = round(timeSect);

% Find the maximum value
maxVal = max(heatValMat(:));

% Initialize an image matrix
heatMap = zeros(1080,1920,3);

% Sets up the color map
colorSize = 128;
mapColor = colormap(jet(colorSize));
close(gcf)

% Loops over all pixels in the heatMap image variable
for k = 1:1920
    for j = 1:1080
        
        % Determines the color
        rgbIdx = floor(heatValMat(j+sizeGauss+1,k+sizeGauss+1)/maxVal *colorSize);
        
        % The minumum index is 1
        if rgbIdx <= 0
            rgbIdx = 1;
        end % END IF
        
        heatMap(j,k,:) = mapColor(rgbIdx,:);
        
    end
end

heatMap(:,:,1) = flipud(heatMap(:,:,1));
heatMap(:,:,2) = flipud(heatMap(:,:,2));
heatMap(:,:,3) = flipud(heatMap(:,:,3));

% plot([640, 640], [1080, 360], 'LineWidth', 10)
% plot([640*2, 640*2], [1080, 360], 'LineWidth', 10)
% plot([640, 640], [200, 0], 'LineWidth', 10)
% plot([640*2, 640*2], [200, 0], 'LineWidth', 10)

heatMap([1:1080-360], [640-25:640+25], :) = zeros(720,51,3);
heatMap([880:1080], [640-25:640+25], :) = zeros(201,51,3);

heatMap([1:1080-360], [640*2-25:640*2+25], :) = zeros(720,51,3);
heatMap([880:1080], [640*2-25:640*2+25], :) = zeros(201,51,3);

heatMapSmall = imresize(heatMap, 0.25);
figure(10)
imshow(heatMapSmall)
imwrite(heatMapSmall,'Laboras_HeatmapExample.png')
end %END FUNCTION

% EOF