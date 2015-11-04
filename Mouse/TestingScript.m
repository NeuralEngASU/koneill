

% img13.time=[]; img13.data=[]; img13.control=[];
% img13.time=img13.time'; img13.control=img13.control';img13.data=permute(img13.data, [2,1]);
% save('img13', 'img13');

for i=1:7
    img7.data(i,:) = img7.data(i,:)./img7.control;
end

for i=1:4
    img9.data(i,:) = img9.data(i,:)./img9.control;
end

for i=1:3
    img13.data(i,:) = img13.data(i,:)./img13.control;
end


%% img7

for i = 1:7
    
    figure(i)
%     plot(1:length(img7.data(1,:)), img7.data(i,:))
    plot(img7.data(i,:))

%     plot(abs(img7.data(i, 1:106)./img7.control(1:106)))
    
end

%% img9

for i = 1:4
    
    figure(i)
%     plot(1:length(img7.data(1,:)), img7.data(i,:))
    plot(img9.data(i,:))

%     plot(abs(img7.data(i, 1:106)./img7.control(1:106)))
    
end


%% img13


for i = 1:3
    
%     figure(i)
%     plot(1:length(img7.data(1,:)), img7.data(i,:))
%     plot(img13.data(i,:))

%     plot(abs(img7.data(i, 1:106)./img7.control(1:106)))

end

% meanData = mean(img13.data, 1);
% 
% plot(meanData);
% 







for i=1:3
    img13.data(i,:) = img13.data(i,:)./img13.control;
end










