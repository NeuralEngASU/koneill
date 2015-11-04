%% Calc postimplant day
ImplantDate = '20130322';
d = regexp(PCA.NEVFile,'\\(\d+)-\d+\\','tokens'); d = cell2mat([d{:}]);
pid = etime(datevec(d,'yyyymmdd'),datevec(ImplantDate,'yyyymmdd'))/60/60/24;


%% PCA classification results nchoosek boxplot
fh = figure('position',[50,50,900,700]);
ah = subplot(1,1,1,'parent',fh);
n = size(PCA.PTotals,2);
boxplot(ah,PCA.PTotals)
ylim(ah,[0,1.02])
hold on
plot(ah,[(1:n)-0.2;(1:n)+0.2],repmat(1./(2:n+1),2,1),'g')
hold off
set(ah,'xtick',1:n,'xticklabel',2:n+1,'ytick',0:0.1:1,'yticklabel',0:10:100)
xlabel('# of classes')
ylabel('Percent correctly classified')
title(sprintf('PCA-LDA classification of %0.0f finger movements and rest (day %0.0f)\n(Electrodes %s)',n,pid,regexprep(num2str(PCA.Electrodes'),'\s+',',')))


%% Plot pairwise confusion matrix
a = nchoosek(PCA.Results.ClassIdxs,2);
b = fliplr(a);
a = a(:,1)+(a(:,2)-1)*length(PCA.Results.ClassIdxs);
b = b(:,1)+(b(:,2)-1)*length(PCA.Results.ClassIdxs);
c = zeros(length(PCA.Results.ClassIdxs)); c(a) = PCA.PTotals(1:size(a,1),1); c(b) = PCA.PTotals(1:size(a,1),1);
n = size(PCA.PTotals,2);

figure('position',[50,50,900,700]);
subplot(1,1,1)
imagesc(c)
cbH = colorbar;
cbPos = get(get(cbH,'ylabel'),'position'); cbPos(1) = cbPos(1)+1;
set(get(cbH,'ylabel'),'string','Percent correctly classified','rotation',270,'position',cbPos)
set(gca,'xtick',1:n+1,'xticklabel',[],'ytick',1:n+1,'yticklabel',PCA.ClassList)
for k=1:n+1
    text(k,n+1.75,PCA.ClassList{k},'horizontalalignment','right','rotation',45)
end
axis square
title(sprintf('Pairwise PCA-LDA classification results\n(Electrodes %s)',regexprep(num2str(PCA.Electrodes'),'\s+',',')))


%% 1 of (9/14) PCA results confusion matrix
n = size(PCA.PTotals,2);
R = PCA.Results.TestResults;
ID = PCA.Results.TestMatID(:,1);
cIdxs = PCA.Results.ClassIdxs;
RMat = nan(10,length(cIdxs));
for k=1:length(cIdxs)
    r = R(ID==cIdxs(k));
    RMat(1:length(r),k) = r;
end

a = permute(repmat(ones(10,1)*(1:n+1),[1,1,n+1]),[1,3,2]);
b = repmat(RMat,[1,1,n+1]);
c = squeeze(sum(a==b,1))./repmat(sum(squeeze(sum(a==b,1)),2),1,n+1);

figure('position',[50,50,900,700]);
subplot(1,1,1)
imagesc(c)
cbH = colorbar;
cbPos = get(get(cbH,'ylabel'),'position'); cbPos(1) = cbPos(1)+1;
set(get(cbH,'ylabel'),'string','Percent classified','rotation',270,'position',cbPos)
set(gca,'xtick',1:n+1,'xticklabel',[],'ytick',1:n+1,'yticklabel',PCA.ClassList)
for k=1:n+1
    text(k,n+1.75,PCA.ClassList{k},'horizontalalignment','right','rotation',45)
end
axis square
ylabel('Classes'); xH = xlabel('Classification Results'); xPos = get(xH,'position'); xPos(2) = xPos(2)+2; set(xH,'position',xPos)
title(sprintf('PCA-LDA 1 of %0.0f classification results (%0.1f%% overall correct)\n(Electrodes %s)',n+1,mean(diag(c))*100,regexprep(num2str(PCA.Electrodes'),'\s+',',')))


%% Plot raw pca results
figure
ClassColors = lines(length(PCA.ClassList));

TrainMat = PCA.Results.TrainMat;
TrainMatID = PCA.Results.TrainMatID;
TestMat = PCA.Results.TestMat;
TestMatID = PCA.Results.TestMatID;
PCTestMat = PCA.Results.PCTestMat;
PCTrainMat = PCA.Results.PCTrainMat;
ClassIdxs = PCA.Results.ClassIdxs;
ClassList = PCA.ClassList;
TrialSamples = size(PCA.Features,4);
L_Results = PCA.Results.TestResults==PCA.Results.TestMatID(:,1);
L_Total = PCA.PTotals(1,end);
e = PCA.Electrodes;

ax(1) = subplot(2,3,[1,4]);
PlotData = (TrainMat - repmat(linspace(0,size(TrainMat,1),size(TrainMat,1))',1,size(TrainMat,2)));
hold on
for k=1:size(TrainMatID,1)
    plot(PlotData(k,:),'color',ClassColors(TrainMatID(k,1)-(min(ClassIdxs)-1),:))
end
plot(repmat((1:length(e)-1)*TrialSamples,2,1),repmat([max(PlotData(:));min(PlotData(:))],1,length(e)-1),'--k')
hold off
axis tight
set(ax(1),'xtick',TrialSamples/2:TrialSamples:TrialSamples*length(e),'xticklabel',e,'ytick',fliplr(-5:-10:-length(ClassIdxs)*10),'yticklabel',fliplr(ClassList(ClassIdxs)),'box','on')
xlabel(ax(1),'Electrode'); ylabel(ax(1),'Movement Features (Rate)'); title(ax(1),'Training')

ax(2) = subplot(2,3,[2,5]);
PlotData = (TestMat - repmat(linspace(0,size(TestMat,1),size(TestMat,1))',1,size(TestMat,2)));
hold on
for k=1:size(TestMatID,1)
    if L_Results(k)
        plot(PlotData(k,:),'color',ClassColors(TestMatID(k,1)-(min(ClassIdxs)-1),:))
    else
        plot(PlotData(k,:),':','color',ClassColors(TestMatID(k,1)-(min(ClassIdxs)-1),:))
    end
end
plot(repmat((1:length(e)-1)*TrialSamples,2,1),repmat([max(PlotData(:));min(PlotData(:))],1,length(e)-1),'--k')
hold off
axis tight
set(ax(2),'xtick',TrialSamples/2:TrialSamples:TrialSamples*length(e),'xticklabel',e,'ytick',[],'yticklabel',[],'box','on')
xlabel(ax(2),'Electrode'); title(ax(2),['Testing (',num2str(L_Total,'%0.2f'),')'])
linkaxes([ax(1),ax(2)],'xy')

ax(3) = subplot(2,3,3);
hold on
for k=1:size(TrainMatID,1)
    plot(PCTrainMat(k,1),PCTrainMat(k,2),'.','color',ClassColors(TrainMatID(k,1)-(min(ClassIdxs)-1),:));
end
hold off
axis tight
xlabel(ax(3),'PC1'); ylabel(ax(3),'PC2'); title(ax(3),'Training')
set(ax(3),'box','on')

ax(4) = subplot(2,3,6);
hold on
for k=1:size(TestMatID,1)
    if L_Results(k)
        plot(PCTestMat(k,1),PCTestMat(k,2),'.','color',ClassColors(TestMatID(k,1)-(min(ClassIdxs)-1),:));
    else
        plot(PCTestMat(k,1),PCTestMat(k,2),'*','color',ClassColors(TestMatID(k,1)-(min(ClassIdxs)-1),:));
    end
end
hold off
axis tight
xlabel(ax(4),'PC1'); ylabel(ax(4),'PC2'); title(ax(4),['Testing (',num2str(L_Total,'%0.2f'),')'])
set(ax(4),'box','on')

xlm = [min([get(ax(3),'xlim'),get(ax(4),'xlim')]),max([get(ax(3),'xlim'),get(ax(4),'xlim')])];
ylm = [min([get(ax(3),'ylim'),get(ax(4),'ylim')]),max([get(ax(3),'ylim'),get(ax(4),'ylim')])];
set(ax(3),'xlim',xlm,'ylim',ylm); set(ax(4),'xlim',xlm,'ylim',ylm)


%% Electrode responses to all classes boxplot

Features = PCA.Features;
Baselines = PCA.Baselines;

% e = findDrivenElects(Features,Baselines);
e = [36,84,12,41];
eList = {'Thumb Flexion','Index Extension','General Flexion/Abduction','General Extension/Abduction'};
ClassList = {'ThumbFlex','IndexFlex','MiddleFlex','RingFlex','LittleFlex','IndexAbd','RingAbd','LittleAbd','ThumbExt','IndexExt','MiddleExt','RingExt','LittleExt','Rest'};

figure('position',[50,50,900,700]);
for m=1:length(e)
    a = [];
    aH = subplot(2,2,m);
    for k=1:length(ClassList)-1
        a(:,k) = squeeze(max(Features(k,e2c(e(m),'pns'),:,:),[],4))*50;
    end
    a(:,k+1) = squeeze(max(Baselines(k,e2c(e(m),'pns'),:,:),[],4))*50;
    bH = boxplot(squeeze(a),'parent',aH);
    title(aH,num2str(e(m)))
    set(aH,'xtick',1:14,'xticklabel',repmat({''},1,14),'ylim',[0-max(a(:))*0.05,max(a(:))+max(a(:))*0.05])
    for k=1:14
        text(k,-max(a(:))/12,ClassList{k},'horizontalalignment','right','rotation',45,'parent',aH)
    end
    set(get(aH,'ylabel'),'string','Rate (Hz)')    
    title(aH,sprintf('%s\n(Electrode %0.0f)',eList{m},e(m)))
end


%% PCA over time

% NEVList{1} = 'D:\Tyler\data\PNS\P201301\20130325-112324\PCA_2000ms_SortedDElects_20130325-112324-001_Sorted_TS.mat';
% NEVList{2} = 'D:\Tyler\data\PNS\P201301\20130329-122625\PCA_2000ms_UnsortedDElects_20130329-122625-001.mat'; %is combined with 20130329-130205
% NEVList{3} = 'D:\Tyler\data\PNS\P201301\20130404-114814\PCA_2000ms_UnsortedDElects_20130404-114814-001.mat';
% NEVList{4} = 'D:\Tyler\data\PNS\P201301\20130408-111310\PCA_2000ms_SortedDElects_20130408-111310-001_Sorted_TS.mat';
% NEVList{5} = 'D:\Tyler\data\PNS\P201301\20130411-112540\PCA_2000ms_SortedDElects_20130411-112540-001_Sorted_TS.mat';
% NEVList{6} = 'D:\Tyler\data\PNS\P201301\20130415-120054\PCA_2000ms_SortedDElects_20130415-120054-001_Sorted_TS.mat';
% NEVList{7} = 'D:\Tyler\data\PNS\P201301\20130417-141736\PCA_2000ms_SortedDElects_20130417-141736-001_Sorted_TS.mat';

NEVList{1} = 'D:\Tyler\data\PNS\P201301\20130325-112324\PCA_2000ms_SortedDElects_20130325-112324-001.mat';
NEVList{2} = 'D:\Tyler\data\PNS\P201301\20130329-122625\PCA_2000ms_UnsortedDElects_20130329-122625-001.mat'; 
NEVList{3} = 'D:\Tyler\data\PNS\P201301\20130404-114814\PCA_2000ms_UnsortedDElects_20130404-114814-001.mat';
NEVList{4} = 'D:\Tyler\data\PNS\P201301\20130408-111310\PCA_2000ms_SortedDElects_20130408-111310-001.mat';
NEVList{5} = 'D:\Tyler\data\PNS\P201301\20130411-112540\PCA_2000ms_SortedDElects_20130411-112540-001.mat';
NEVList{6} = 'D:\Tyler\data\PNS\P201301\20130415-120054\PCA_2000ms_SortedDElects_20130415-120054-001.mat';
NEVList{7} = 'D:\Tyler\data\PNS\P201301\20130417-141736\PCA_2000ms_SortedDElects_20130417-141736-001.mat';

ImplantDate = '20130322';

RMat9 = nan(10000,length(NEVList));
RMat14 = nan(10000,length(NEVList));
for m=1:length(NEVList)
    d = regexp(NEVList{m},'\\(\d+)-\d+\\','tokens'); d = cell2mat([d{:}]);
    pid(m) = etime(datevec(d,'yyyymmdd'),datevec(ImplantDate,'yyyymmdd'))/60/60/24;
    load(NEVList{m})    
    if m<=2
%         r = PCA.Results.TestResults==PCA.Results.TestMatID(:,1);
%         R9 = zeros(10000,1);
%         for k=1:10000
%             rIdx = randi(length(r),[length(r),1]);
%             R9(k) = sum(r(rIdx))/length(r);
%         end
        R9 = PCA.PTotals(:,8);
        R14 = nan;
    else
        R9 = PCA.PTotals(:,8);
        R14 = PCA.PTotals(:,13);
%         r = PCA.Results.TestResults==PCA.Results.TestMatID(:,1);
%         R14 = zeros(10000,1);
%         for k=1:10000
%             rIdx = randi(length(r),[length(r),1]);
%             R14(k) = sum(r(rIdx))/length(r);
%         end
    end    
    RMat9(1:length(R9),m) = R9;    
    RMat14(1:length(R14),m) = R14;    
end

%%
fH = figure('position',[50,50,1000,600]); 
aH = subplot(1,1,1,'parent',fH);
set(aH,'nextplot','add')
boxplot(RMat9,'widths',0.2,'symbol','','parent',aH)
boxplot(RMat14,'widths',0.2,'symbol','','parent',aH)
ylim(aH,[0,1.1])
set(aH,'xtick',1:length(pid),'xticklabel',pid)
xlabel('Postimplant Day')
ylabel('Percent correctly classified')
title('PCA-LDA Classification Over Time')

hg = get(aH,'children');
c1 = get(hg(1),'children');
for k=1:length(c1) %14DOF
    if strcmp(get(c1(k),'tag'),'Median')
        set(get(get(c1(k),'Annotation'),'LegendInformation'),'IconDisplayStyle','on')
    else
        set(get(get(c1(k),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    end
    if ~isempty(get(c1(k),'tag'))
        set(c1(k),'xdata',get(c1(k),'xdata')+0.15,'color','k') 
    else
        delete(c1(k));
    end
end

c2 = get(hg(2),'children');
[~,na] = unique(get(c2,'tag'));
for k=1:length(c2) %9DOF
    if strcmp(get(c2(k),'tag'),'Median')
        set(get(get(c2(k),'Annotation'),'LegendInformation'),'IconDisplayStyle','on')
    else
        set(get(get(c2(k),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    end
    if ~isempty(get(c2(k),'tag'))
        set(c2(k),'color','b')
        if ~any(k==[na;na-1])
            set(c2(k),'xdata',get(c2(k),'xdata')-0.15)
        end
    else
        delete(c2(k));
    end
end

hA = get(hg,'Annotation');
hLL = get([hA{:}],'LegendInformation');
set([hLL{:}],{'IconDisplayStyle'},{'on','on'}')
set(hg,{'DisplayName'},{'14 Classes','9 Classes'}')
legend('show','location','northwest');