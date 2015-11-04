function [sig,p,stats] = check_pli(grid,srcdir,alpha)
if nargin<3||isempty(alpha),alpha=0.05;end
if nargin<2||isempty(srcdir),srcdir='E:\data\ecogres';end

% load the data
tmp = load(fullfile(srcdir,sprintf('g%dmpli_UNR.mat',grid)));
npair = size(tmp.p,2);

% run first pair to get the output information
tmpp_real = squeeze(tmp.p(:,1,1));
tmpp_surr = squeeze(tmp.p(:,1,2:end));
[~,~,st] = kruskalwallis(...
    [tmpp_surr(:); tmpp_real(:)],... % all data concatenated
    [zeros(numel(tmpp_surr),1); ones(numel(tmpp_real),1)],... % group identifiers
    'off'); % no display
[c,m] = multcompare(st,'Alpha',alpha,'Display','off');

sig = nan(1,npair); sig(1)=sign(m(2,1)-m(1,1))*(c(end)<alpha);
p = nan(1,npair); p(1)=c(end);
stats = nan(npair,4); stats(1,:) = [mean(tmpp_real(:)) std(tmpp_real(:)) mean(tmpp_surr(:)) std(tmpp_surr(:))];
for kk=2:npair
    tmpp_real = squeeze(tmp.p(:,kk,1));
    tmpp_surr = squeeze(tmp.p(:,kk,2:end));
    [~,~,st] = kruskalwallis(...
        [tmpp_surr(:); tmpp_real(:)],... % all data concatenated
        [zeros(numel(tmpp_surr),1); ones(numel(tmpp_real),1)],... % group identifiers
        'off'); % no display
    [c,m] = multcompare(st,'Alpha',alpha,'Display','off');
    sig(kk)=sign(m(2,1)-m(1,1))*(c(end)<alpha);
    p(kk)=c(end);
    stats(kk,:) = [mean(tmpp_real(:)) std(tmpp_real(:)) mean(tmpp_surr(:)) std(tmpp_surr(:))];
end

%% plots
figw = min(max(300,90+0.4*npair),1850);
figure('Position',[50 200 figw 600],'PaperPositionMode','auto','Name',sprintf('Grid %d PLI',grid),'NumberTitle','off');
ax(1) = axes('Units','pixels','Position',[60 335 figw-(50+25) 230]);
order = [1 3];
if sum(stats(:,3))>sum(stats(:,1)),order=[3 1];end
b = bar(stats(:,order),'EdgeColor','none');
if order(1)==1
    b(1).FaceColor = [0 0.447 0.741];
    b(1).EdgeColor = [0 0.447 0.741];
    b(2).FaceColor = [0.85 0.325 0.098];
    b(2).EdgeColor = [0.85 0.325 0.098];
else
    b(2).FaceColor = [0 0.447 0.741];
    b(2).EdgeColor = [0 0.447 0.741];
    b(1).FaceColor = [0.85 0.325 0.098];
    b(1).EdgeColor = [0.85 0.325 0.098];
end
xlim([1 npair]);
title('Group Means');
ylabel('PLI');
if order(1)==1
    legend({'Real','Surrogate'});
else
    legend({'Surrogate','Real'});
end
ax(2) = axes('Units','pixels','Position',[60  50 figw-(50+25) 230]);
b = bar(sig,'EdgeColor','none');
b.FaceColor = [0.3 0.3 0.3];
b.EdgeColor = [0.3 0.3 0.3];
xlim([1 npair]);
ylim([-1.1 1.1]);
title('Significance');
xlabel('Channel Pair');
ylabel('PLI');
linkaxes(ax,'x');