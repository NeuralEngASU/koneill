% generate figures of average cross-correlation function vs. separation distance
clear all;
close all;
clc;

% set up environment
% use('cleanpath');
% use('pkgs/skellis');

% define which grids to pull in
grids=defgrids;

% specify rereferencing
ref='unr';

% loop over the grids
for g=1:length(grids)
    load(['D:\Kevin\SpencerPaper\201203\g' num2str(g) 'mccf_' upper(ref) '.mat'],'mccf','chanpairs','lags','spacing','layout');

    % get average
    [mccf,seps,pairsep]=avgpersep(mccf,chanpairs,layout,spacing,grids(g).badchan);

    % only use a subset of lags to create the mesh
    lagidx=fix(logspace(0,log10((length(lags)-1)/2),30));
    lagidx=lagidx+(length(lags)-1)/2;
    lagidx=sort(unique(lagidx),'ascend');
    sidx=fix(linspace(1,length(seps),20));
    sidx=sort(unique(sidx),'ascend');

    % open figure
    figure('Name',['Cross-correlation: ' grids(g).label],'PaperPositionMode','auto','Position',[140 237 1676 375]);
    cmap=colormap('copper');

    % TITLE
    ax(1)=subplot('Position',[0 0.07 0.03 0.89]);
    set(ax(1),'Visible','off');
    text(0.5,0.5,grids(g).shortname,'Parent',ax(1),'FontSize',14,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle','Rotation',90);

    % MESH 3D PLOT OF AVG CCF
    ax(2)=subplot('Position',[0.06 0.13 0.17 0.8]);

    % plot the lines in the "lag" direction
    plot3(ax(2),repmat(lags(lagidx),1,length(seps(sidx))),...
        repmat(seps(sidx)',length(lagidx),1),...
        mccf(lagidx,sidx),'b')
    hold on

    % plot the lines in the "sep" direction
    plot3(ax(2),repmat(lags(lagidx),1,length(seps(sidx)))',...
        repmat(seps(sidx)',length(lagidx),1)',...
        mccf(lagidx,sidx)','b');
    hold off;

    % set display properties
    xlim([floor(lags(lagidx(1))) ceil(lags(lagidx(end)))]);
    ylim([seps(sidx(1))-0.2 seps(sidx(end))+0.2]);
    zlim([-0.2 1.1]);
    zlabel('correlation');
    set(ax(2),'XGrid','on','YGrid','on','ZGrid','on');
    set(ax(2),'XScale','log','YScale','lin','ZScale','lin');
    set(ax(2),'XTick',[logspace(-3,0,4) 2],'XTickLabel',{'0.001','','','','2'});
    yticks=linspace(0,floor(seps(sidx(end))),5);
    yticklabels=cell(size(yticks));
    yticklabels(:)={' '};
    yticklabels{1}=round(seps(1));
    yticklabels{end}=floor(seps(sidx(end)));
    set(ax(2),'YTick',yticks,'YTickLabel',yticklabels)
    set(ax(2),'ZTick',[0 0.5 1],'ZTickLabel',{'0','','1'});
    set(ax(2),'TickLength',[0 0]);
    view([115 49]);

    % PROJECTION TO PLOT CCF VS. LAG, DIFF LINE PER SEP
    ax(3)=subplot('Position',[0.285 0.15 0.28 0.7]);
    set(ax(3),'ColorOrder',cmap(fix(linspace(1,64,length(sidx))),:));
    hold all;
    sep_legstr=cell(1,length(sidx));
    for k=1:length(sidx)
        plot(lags,mccf(:,k),'LineWidth',2);
        sep_legstr{k}=[num2str(seps(sidx(k)),2) ' mm'];
    end
    llegstr=repmat({' '},1,64);
    llegstr(fix(linspace(1,64,length(sidx))))=sep_legstr;
    lcolorbar(llegstr,'Parent',ax(3));
    xlabel('lag (s)');
    xlim([floor(lags(lagidx(1))) ceil(lags(lagidx(end)))]);
    ylim([-0.2 1.1]);
    % set(ax(2),'XScale','log');
    grid on;

    % add inset to highlight dynamics at the beginning of the x-axis
    ax(4)=axes('Position',[0.32 0.4 0.20 0.4]);
    set(ax(4),'ColorOrder',cmap(fix(linspace(1,64,length(sidx))),:));
    hold all;
    for k=1:length(sidx)
        plot(lags,mccf(:,k),'LineWidth',2);
    end
    xlim([floor(lags(lagidx(1))) ceil(lags(lagidx(4)))]);
    ylim([-0.2 1.1]);
    set(ax(4),'XScale','log');
    set(ax(4),'XTick',logspace(-3,0,4),'XTickLabel',{'0.001','','','1'});
    set(ax(4),'YTick',[0 0.5 1],'YTickLabel',{'0','','1'});
    set(ax(4),'TickLength',[0 0]);
    box on;
    grid on;

    % PROJECT TO PLOT CCF VS. SEP, DIFF LINE PER LAG
    ax(5)=subplot('Position',[0.65 0.15 0.3 0.7]);
    set(ax(5),'ColorOrder',cmap(fix(linspace(1,64,length(lagidx))),:));
    hold all;
    for k=1:length(lagidx)
        plot(seps,mccf(lagidx(k),:),'LineWidth',2);
    end
    idx=1;
    lag_legstr=cell(1,round(length(lagidx)/3));
    for k=fix(linspace(1,length(lagidx),10))
        lag_legstr{idx}=[num2str(lags(lagidx(k)),2) ' sec'];
        idx=idx+1;
    end
    llegstr=repmat({' '},1,64);
    llegstr(fix(linspace(1,64,length(lag_legstr))))=lag_legstr;
    lcolorbar(llegstr,'Parent',ax(5));
    xlabel('separation (mm)');
    xlim([seps(sidx(1)) seps(sidx(end))]);
    ylim([-0.2 1.1]);
    grid on;

    saveas(gcf,['D:\Kevin\SpencerPaper\201203\figures\' grids(g).shortname '_' upper(ref) '.fig']);
    saveas(gcf,['D:\Kevin\SpencerPaper\201203\figures\' grids(g).shortname '_' upper(ref) '.png']);
    saveas(gcf,['D:\Kevin\SpencerPaper\201203\figures\' grids(g).shortname '_' upper(ref) '.eps'],'epsc2');

%     close(gcf);
end