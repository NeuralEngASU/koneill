load('DeltaSpeech_Day1_DS5000_BandPass250_2000_WPLI_winSize0.1.mat')
load('DeltaSpeech_Day1_DS5000_BandPass250_2000.mat')

%%

wordCoord = floor(pts_words{2}./6./Header.Fs./0.1);

pliData = [];
test = mean(p,2);

for ii = 1:32%:length(wordCoord)
    pliData(:,ii) = p(wordCoord(ii)-10:(wordCoord(ii)+20),17);
    testData(:,ii) = test(wordCoord(ii)-10:(wordCoord(ii)+20));
end

t = -1:0.1:2;

meanPliData = mean(pliData,2);
meanTestData = mean(testData,2);

plot(t,meanPliData);
ylim([0,1])


%%
figure
for ii = 32:51
    plot([wordCoord(ii), wordCoord(ii)]*0.1, [0,1], 'k')
    hold on
end
    
plot([7667:8019]*0.1, p(7667:8019,1))
xlim([766, 802])
ylim([-0.5,1.5])

% EOF