% compute coherence
clear all;
close all;
clc;

% create parallel pool
if isempty(gcp('nocreate')),parpool(16);end

% data segmenting
movingwin=[2 2];

% number of surrogate data iterations
nsur = 100;

% specify rereferencing
ref='unr';

disp(' ');
disp('GEN_PLI: Generate phase-lag index');
disp(['Rereferencing scheme: ' ref]);
disp(' ');

% load data
grids=defgrids;
srcdir='E:\data\PLI\DownsampledData';%'S:\Spencer\ecogres';
tgtdir='E:\data\PLI\DownsampledData';%'S:\Spencer\ecogres';
for g=1:length(grids)
    
    % update user
    disp(['grid ' num2str(g) '/' num2str(length(grids))]);
    layout=grids(g).layout;
    spacing=grids(g).spacing;
    
    % load data
    disp('  loading data');
    [raw,params]=get_segdata(srcdir,grids(g),g,movingwin,ref);
    nchan=size(raw,3);
    nwin=size(raw,2);
    
    % set up channel pairs
    chanpairs=nchoosek(sort(unique(layout(~isnan(layout)&layout>0)),'ascend'),2);
    npair=size(chanpairs,1);
    
    % initialize variables
    p=zeros(nwin,npair,nsur+1);
    r=zeros(nwin,npair,nsur+1);
    
    % for surrogate data: random phase shift applied to each channel
    smp = [zeros(npair,1) 500+round(500*rand(npair,nsur))]; % 500-1000 samples (0.5 - 1.0 seconds)
    
    % now run the rest of the pairs
    disp('  processing');
    tocs = nan(1,npair);
    for cp=1:npair
        stopwatch = tic;
        remaining = (npair-cp+1)*nanmedian(tocs);
        remstr = 'waiting to estimate time remaining';
        if ~isnan(remaining), remstr = sprintf('%s remaining',hms(remaining,'hh:mm:ss')); end%Utilities.hms(remaining,'hh:mm:ss')); end
        fprintf('\tchanpair %d/%d (%s)\n',cp,size(chanpairs,1),remstr);
        
        % pull out data for this iteration
        raw1=raw(:,:,chanpairs(cp,1));
        raw2=raw(:,:,chanpairs(cp,2));
        
        % calculate PLI and R for real data
        tmpp = nan(nwin,nsur+1);
        tmpr = nan(nwin,nsur+1);
        tmpsmp = smp(cp,:);
        parfor ss=1:nsur+1
            
            % circshift the second channel (destroy correlations)
            raw2s = circshift(raw2,[tmpsmp(ss) 0]);
            
            % calculate PLI and R
            [tmpp(:,ss),tmpr(:,ss)]=pli(raw1, raw2s);%util.pli(raw1,raw2s);
        end
        p(:,cp,:) = tmpp;
        r(:,cp,:) = tmpr;
        tocs(cp) = toc(stopwatch);
    end
    
    % save results
    save(fullfile(tgtdir,sprintf('g%dmpli_%s.mat',g,upper(ref))),'p','r','chanpairs','spacing','layout','smp');
end