function [ imgOut ] = GenHeatmap( data, saveData )

% Normalize Data
minData = min(data(:));
data = data + abs(minData);
maxData = max(data(:));

% Select colormap and aquire colormap values
colormap(jet)
mapColor = jet;
close(gcf)

% Initialization of image
imgOut = zeros(2*size(data,1), size(data,2)/2, 3); 

% Loop over all rows and half of the columns of data
for i=1:size(data,1)
    for j=1:size(data,2)/2
%         fprintf('%d, %d\n', i, j)
        
        % Calculate the rgb index for the current data values
        rgbIdx = floor(mean(data(i,2*j-1:2*j))/maxData *64);
        
        % The minumum index is 1
        if rgbIdx <= 0
            rgbIdx = 1;
        end % END IF
        
        % set two rows of the output image to be the same color
        imgOut(i*2, j, :) = mapColor(rgbIdx,:);
        imgOut(i*2-1,j,:) = imgOut(i*2,j,:);
        
        
    end % END FOR
end % END FOR

imshow(imgOut)

shortName = 'Heatmap';

if saveData
    imwrite(imgOut, ['G:\CodeRepo\MATLAB\Mouse\Other\', shortName, '.png']);
    saveas(gcf,['G:\CodeRepo\MATLAB\Mouse\Other\', shortName, '.fig']);
    saveas(gcf,['G:\CodeRepo\MATLAB\Mouse\Other\', shortName, '.eps'],'epsc2');
end % END IF

end % END FUNCTION

% EOF