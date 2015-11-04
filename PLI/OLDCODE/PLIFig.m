

%% Load Data
% 
% load('PLI.mat')
% % 
% PLI = PLI2;
% R = R2;
% zScore = zScoreReport2;

%% Copy half of confusion matrix

R2 = R;
R3 = permute(R, [2,1,3]);
R2(R2 == 0) = R3(R2==0);

R = R2;

PLI2 = PLI;
PLI3 = permute(PLI, [2,1,3]);
PLI2(PLI2 == 0) = PLI3(PLI2==0);

PLI = PLI2;

zScore2 = zScore;
zScore3 = permute(zScore, [2,1,3]);
zScore2(zScore2 == 0) = zScore3(zScore2==0);

zScore = zScore2;

<<<<<<< HEAD
=======
folder = 'delta';

>>>>>>> origin/master
%% Plot Basic Figures

% Plot R values for real real signals and average for surrogate signals

figure(1)
colormap('jet')
imagesc(R(:,:,1), [0,1]);

set(gca, 'XTick', []);
<<<<<<< HEAD
set(gca, 'XTickLabels', []);
set(gca, 'YTick', []);
set(gca, 'YTickLabels', []);
colorbar(gca)

title('R Matrix, phase difference between signals')

=======
set(gca, 'XTickLabel', []);
set(gca, 'YTick', []);
set(gca, 'YTickLabel', []);
colorbar()

title('R Matrix, phase difference between signals')

print(['Results\', folder, '\RReal'],'-dpng')
savefig(['Results\', folder, '\RReal'])

>>>>>>> origin/master
figure(2)
colormap('jet')
imagesc(mean(R(:,:,2:end),3), [0,1]);

set(gca, 'XTick', []);
<<<<<<< HEAD
set(gca, 'XTickLabels', []);
set(gca, 'YTick', []);
set(gca, 'YTickLabels', []);
colorbar(gca)

title('R Matrix for Surrogate Data')

=======
set(gca, 'XTickLabel', []);
set(gca, 'YTick', []);
set(gca, 'YTickLabel', []);
colorbar()

title('R Matrix for Surrogate Data')

print(['Results\', folder, '\RSurrogate'],'-dpng')
savefig(['Results\', folder, '\RSurrogate'])

>>>>>>> origin/master
%% Plot PLI and average PLI for surrogate data

figure(3)
colormap('jet')
imagesc(PLI(:,:,1), [0,0.25]);

set(gca, 'XTick', []);
<<<<<<< HEAD
set(gca, 'XTickLabels', []);
set(gca, 'YTick', []);
set(gca, 'YTickLabels', []);
colorbar(gca)

title('Phase Lag Index')

=======
set(gca, 'XTickLabel', []);
set(gca, 'YTick', []);
set(gca, 'YTickLabel', []);
colorbar()

title('Phase Lag Index')

print(['Results\', folder, '\PLIReal'],'-dpng')
savefig(['Results\', folder, '\PLIReal'])

>>>>>>> origin/master
figure(4)
colormap('jet')
imagesc(mean(PLI(:,:,2:end),3), [0,0.25]);

set(gca, 'XTick', []);
<<<<<<< HEAD
set(gca, 'XTickLabels', []);
set(gca, 'YTick', []);
set(gca, 'YTickLabels', []);
colorbar(gca)

title('Phase Lag Index, Surrogate Data')

=======
set(gca, 'XTickLabel', []);
set(gca, 'YTick', []);
set(gca, 'YTickLabel', []);
colorbar()

title('Phase Lag Index, Surrogate Data')

print(['Results\', folder, '\PLISurrogate'],'-dpng')
savefig(['Results\', folder, '\PLISurrogate'])

>>>>>>> origin/master
%% zScore for data and surrogate data

figure(5)
colormap('jet')
imagesc(zScore(:,:,1), [-5,5]);

set(gca, 'XTick', []);
<<<<<<< HEAD
set(gca, 'XTickLabels', []);
set(gca, 'YTick', []);
set(gca, 'YTickLabels', []);
colorbar(gca)

title('Z-Score')

=======
set(gca, 'XTickLabel', []);
set(gca, 'YTick', []);
set(gca, 'YTickLabel', []);
colorbar()

title('Z-Score')

print(['Results\', folder, '\ZScoreReal'],'-dpng')
savefig(['Results\', folder, '\ZScoreReal'])

>>>>>>> origin/master
figure(6)
colormap('jet')
imagesc(mean(zScore(:,:,2:end),3), [-5,5]);

set(gca, 'XTick', []);
<<<<<<< HEAD
set(gca, 'XTickLabels', []);
set(gca, 'YTick', []);
set(gca, 'YTickLabels', []);
colorbar(gca)

title('Z-Score, Surrogate Data')

=======
set(gca, 'XTickLabel', []);
set(gca, 'YTick', []);
set(gca, 'YTickLabel', []);
colorbar()

title('Z-Score, Surrogate Data')

print(['Results\', folder, '\ZScoreSurrogate'],'-dpng')
savefig(['Results\', folder, '\ZScoreSurrogate'])

>>>>>>> origin/master
%% zScore pass for signifigance

zScorePass = (abs(zScore(:,:,1)) >= 1.96) .* sign(zScore(:,:,1));

figure(7)
colormap('jet')
imagesc(zScorePass, [-1,1]);

set(gca, 'XTick', []);
<<<<<<< HEAD
set(gca, 'XTickLabels', []);
set(gca, 'YTick', []);
set(gca, 'YTickLabels', []);
colorbar(gca)

title(sprintf('Z-Score >= 1.96 (.95 Confidence)\n+1 = Signifigant Coupling. -1 = Signifigant Non-Coupling'))

=======
set(gca, 'XTickLabel', []);
set(gca, 'YTick', []);
set(gca, 'YTickLabel', []);
colorbar()

title(sprintf('Z-Score >= 1.96 (.95 Confidence)\n+1 = Signifigant Coupling. -1 = Signifigant Non-Coupling'))
print(['Results\', folder, '\ZScoreSignificance'],'-dpng')
savefig(['Results\', folder, '\ZScoreSignificance'])


%% PLI Distribution

figure(8)
subplot(2,2,1)
hist(squeeze(PLI(1,1,:))); title(sprintf('C1 vs C1')); xlabel('PLI'); ylabel('Count');
% hist(squeeze(PLI(1,1,:))); title(sprintf('C1 vs C1')); xlabel('PLI'); ylabel('Count');
xlim([0,0.3])

subplot(2,2,2)
hist(squeeze(PLI(13,3,:))); title(sprintf('C13 vs C3, Positive-Significant')); xlabel('PLI'); ylabel('Count');
% hist(squeeze(PLI(12,3,:))); title(sprintf('C12 vs C3, Positive-Significant')); xlabel('PLI'); ylabel('Count');
xlim([0,0.3])

subplot(2,2,3)
hist(squeeze(PLI(23,6,:))); title(sprintf('C23 vs C6, Non-Significant')); xlabel('PLI'); ylabel('Count');
% hist(squeeze(PLI(22,6,:))); title(sprintf('C22 vs C6, Non-Significant')); xlabel('PLI'); ylabel('Count');
xlim([0,0.3])

subplot(2,2,4)
hist(squeeze(PLI(20,7,:))); title(sprintf('C20 vs C7, Negative-Significant')); xlabel('PLI'); ylabel('Count');
% hist(squeeze(PLI(24,5,:))); title(sprintf('C24 vs C5, Negative-Significant')); xlabel('PLI'); ylabel('Count');
xlim([0,0.3])

%% PLI Differrence Distribution

figure(9)
subplot(2,2,1)
hist(squeeze(PLI(1,1,:)-PLI(1,1,1))); title(sprintf('C1 vs C1')); xlabel('PLI'); ylabel('Count');
% hist(squeeze(PLI(1,1,:))); title(sprintf('C1 vs C1')); xlabel('PLI'); ylabel('Count');
% xlim([0,0.3])

subplot(2,2,2)
hist(squeeze(PLI(13,3,:) - PLI(13,3,1))); title(sprintf('C13 vs C3, Positive-Significant')); xlabel('PLI'); ylabel('Count');
% hist(squeeze(PLI(12,3,:))); title(sprintf('C12 vs C3, Positive-Significant')); xlabel('PLI'); ylabel('Count');
% xlim([0,0.3])

subplot(2,2,3)
hist(squeeze(PLI(23,6,:) - PLI(23,6,1))); title(sprintf('C23 vs C6, Non-Significant')); xlabel('PLI'); ylabel('Count');
% hist(squeeze(PLI(22,6,:))); title(sprintf('C22 vs C6, Non-Significant')); xlabel('PLI'); ylabel('Count');
% xlim([0,0.3])

subplot(2,2,4)
hist(squeeze(PLI(20,7,:) - PLI(20,7,1))); title(sprintf('C20 vs C7, Negative-Significant')); xlabel('PLI'); ylabel('Count');
% hist(squeeze(PLI(24,5,:))); title(sprintf('C24 vs C5, Negative-Significant')); xlabel('PLI'); ylabel('Count');
% xlim([0,0.3])
>>>>>>> origin/master
