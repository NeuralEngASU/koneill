% Generates artifical Ca2+ transients

function [ transData, time ] = GenCa2Trans(  )

transData(1,:) = ones(1,500);
transData(2,:) = ones(1,500);
transData(3,:) = ones(1,500);
transData(4,:) = ones(1,500);
transData(5,:) = ones(1,500);

time = 0:(1/100):4.99;
% 
% tmpTime = 0:1:249;
% a = 1./exp(-tmpTime/50);
% plot(1./(exp(-tmpTime/20))/2)
% 
% 
% plot(70*(exp(-tmpTime/50)))

for i = 1:250
    
    randNoise = rand(1)*1.5;
    
    tmpTime = 0:1:99;
    transData(i,51:150) = 1./(exp(-tmpTime./(20+randNoise)))/2;
    
    randNoise = rand(1)*1.5;
    tmpTime = 0:1:299;
    tmpData = transData(i,150) * (exp(-tmpTime/(70+randNoise)));
    transData(i,151:450) = tmpData;
    
end % END FOR


for i = 1:size(transData, 1)
    transData(i,:) = awgn(transData(i,:), 2.5);
end % END FOR

% figure(1)
% hold on
% plot(time,transData(1,:),'r')
% plot(time,transData(2,:),'g')
% plot(time,transData(3,:),'b')
% plot(time,transData(4,:),'k')
% plot(time,transData(5,:),'m')
% hold off

title(sprintf('Raw Fake Ca^{2+} Transient Data\n2.5 snr, tau=20,70'))
xlabel('Time, sec')
ylabel('\DeltaF/F')

set(gcf, 'Units', 'inches')
set(gcf, 'Position',[0 0 5 3])
set(gcf, 'PaperUnits','inches','PaperPosition',[2 2 5 3])


end % END FUNC

% EOF



