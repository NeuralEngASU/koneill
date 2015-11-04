%%
mouseData = round(mouseData);

% sigma = 3; mean = 0;
% x=-10:10;
% fx=1/sqrt(2*pi)/sigma*exp(-(x-mean).^2/2/sigma/sigma);
% 
% plot(x, fx)


%%

sizeGauss = 100;

x=linspace(-sizeGauss, sizeGauss, 2*sizeGauss+1);
y=x';               
[X,Y]=meshgrid(x,y);
z=exp(-(X.^2+Y.^2)/(2*1000));
surf(x,y,z);shading interp




%%

heatValMat = zeros(1080 + 2*sizeGauss+1, 1920 + 2*sizeGauss + 1);  


for i = 1:size(mouseData, 1)
    xPos = mouseData(i,1) + sizeGauss+1;
    yPos = mouseData(i,2) + sizeGauss+1;
    
    yMat = [yPos-sizeGauss:yPos+sizeGauss];
    xMat = [xPos-sizeGauss:xPos+sizeGauss];
    
    heatValMat(yMat, xMat) = heatValMat(yMat, xMat) + z;
    
    
end % END FOR

%%
maxVal = max(heatValMat(:));

heatMap = zeros(1080,1920,3);

colorSize = 128;
mapColor = colormap(hot(colorSize));
close(gcf)

for k = 1:1920
    for j = 1:1080
        
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