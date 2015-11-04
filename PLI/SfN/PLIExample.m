%% PLI Example

t = 0:2*pi/1000:(2*pi-2*pi/1000);

s = sin(t)';

s2 = [];
for ii = 1:1000
    
    s2(:,ii) = circshift(s,ii-1);
    
end % END FOR


s3 = repmat(s2,10,1);
ss = repmat(s,10,1000);

noise1 = randn(10000,1000)./2.5;
noise2 = randn(10000,1000)./2.5;

s3 = s3+noise1;
ss = ss+noise2;

% plot(t,sin(t))
%%
figure(1)
subplot(2,1,1)

plot(ss(1:1000,1))
box on
ylim([-2.5,2.5]);

title('Reference Signal')
xlim([1,1000])
set(gca, 'XTick', [])
set(gca, 'YTick', [])

subplot(2,1,2)

plot(s3(1:1000,1), 'color',[0.8500    0.3250    0.0980])
box on
ylim([-2.5,2.5]);

title('Paired Signal')
xlim([1,1000])
set(gca, 'XTick', [])
set(gca, 'YTick', [])


pos = [7.0729    5.7917    2.0313    4.3854];

set(gcf, 'Units', 'Inches')
set(gcf, 'position', pos)
set(gcf, 'paperposition', pos)

%%

figure(3)
subplot(2,1,1)

plot(circshift(s(1:1000,1),250), 'linewidth', 3,'color',[0.8500    0.3250    0.0980])
box on
ylim([-2.5,2.5]);

title('Reference Signal')
xlim([1,1000])
set(gca, 'XTick', [])
set(gca, 'YTick', [])

subplot(2,1,2)

plot(circshift(s(1:1000,1),500), 'color',[0.8500    0.3250    0.0980],'linewidth', 3)
box on
ylim([-2.5,2.5]);

title('Paired Signal')
xlim([1,1000])
set(gca, 'XTick', [])
set(gca, 'YTick', [])


pos = [7.0729    5.7917    2.0313    4.3854];

set(gcf, 'Units', 'Inches')
set(gcf, 'position', pos)
set(gcf, 'paperposition', pos)

%% Find Instantaneous Phase

sPhase = atan2(imag(hilbert(s3)),s3);

sOrig = atan2(imag(hilbert(ss)),ss);

deltaPhi = sOrig - sPhase;
deltaPhi(deltaPhi < -pi) = deltaPhi(deltaPhi < -pi) + 2*pi;
deltaPhi(deltaPhi >  pi) = deltaPhi(deltaPhi >  pi) - 2*pi;

deltaPhi2 = deltaPhi;
% 
% deltaPhi2(deltaPhi2 < -pi/2) = deltaPhi2(deltaPhi2 < -pi/2) + pi/2;
% deltaPhi2(deltaPhi2 >  pi/2) = deltaPhi2(deltaPhi2 >  pi/2) - pi/2;

p = abs(mean(abs(deltaPhi).*sign(deltaPhi),1))'./mean(abs(deltaPhi),1)';
% p = abs(mean(sign(deltaPhi)));

pPlot = circshift(p,floor(length(p)/2));
phasePlot =linspace(-pi,pi,length(p));

figure(3)

plot(phasePlot,pPlot)
% ylim([0,1])
xlim([-pi,pi+eps])

set(gca, 'XTick', [-pi, -pi/2, 0, pi/2, pi])
str={'-pi';'-pi/2';'0';'pi/2';'pi'};
set(gca,'xticklabel',str)

xlabel('Phase')
ylabel('WPLI')


