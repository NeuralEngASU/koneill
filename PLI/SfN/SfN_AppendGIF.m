%% Append Gif

expr = '(2014PP04Sz4_15Sec_Frame_[0-9]+\.png)';

imFileNames = dir(['D:\PLI\SeizureDetection\AnimatedFigures', '\*.png']);

for ii = 1:40
    listNum = 61:100;
    
    imIdx = listNum(ii);
    
    imData = imread(imFileNames(imIdx).name);
    
    [A,map] = rgb2ind(imData,256);
    if ii == 1;
        imwrite(A,map,'2014PP04Sz4_15Sec_gif.gif','gif','LoopCount',Inf,'DelayTime',0.5);
    else
        imwrite(A,map,'2014PP04Sz4_15Sec_gif.gif','gif','WriteMode','append','DelayTime',0.5);
    end
    
end % END FOR