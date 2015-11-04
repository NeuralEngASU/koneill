% generate data
clear all;
close all;
clc;

% set up environment

% establish data sources
sources(1).marker='E:\data\human CNS\Utah\delta\verbal\MARKERS_day1.mat';
sources(1).nsfile='E:\data\human CNS\Utah\delta\verbal\20080730-151259\20080730-151259-002.ns5';
sources(2).marker='E:\data\human CNS\Utah\delta\verbal\MARKERS_day2.mat';
sources(2).nsfile='E:\data\human CNS\Utah\delta\verbal\20080731-111811\20080731-111811-001.ns5';
sources(3).marker='E:\data\human CNS\Utah\delta\verbal\MARKERS_day3.mat';
sources(3).nsfile='E:\data\human CNS\Utah\delta\verbal\20080801-102946\20080801-102946-001.ns5';
sources(4).marker='E:\data\human CNS\Utah\delta\verbal\MARKERS_day4.mat';
sources(4).nsfile='E:\data\human CNS\Utah\delta\verbal\20080801-135921\20080801-135921-001.ns5';

fullwin = [-1.5, 2];
newFs = 5000; % Downsample to this

% % some filter characteristics
% Apass=1;   % Passband Ripple (dB)
% Astop=20;  % Stopband Attenuation (dB)

% loop over sources
data  = cell(1,length(sources));
class = cell(1,length(sources));
for s = 1:length(sources)

    % update user
    disp(['Source ' num2str(s) '/' num2str(length(sources))]);

    % open NS5 and marker data
    Header = openNSx(sources(s).nsfile);
    Header = Header.MetaTags;
    Fs = Header.SamplingFreq;
    ns.fs = Fs;
   
    load(sources(s).marker,'offset_words','pts_words','lbl_words');

    % set up indices to read data from
    tmpclass = [];
    for c = 1:length(pts_words)
        tmpclass = cat(2, tmpclass, repmat(c, 1, length(pts_words{c})));
    end
    [allpts,sidx] = sort(cat(2,pts_words{:}));
    tmpclass      = tmpclass(sidx);
    offsets       = cat(2,offset_words{:});
    offsets       = offsets(sidx);
    starts        = allpts + fix(fullwin(1)*Fs);
    total         = fix(diff(fullwin)*Fs);

    % loop through trial markers
    data{s} = zeros((total)/(Fs/newFs), 32, length(allpts));
    
    for p=1:length(allpts)

        % read extra to discard for filter start-up effects
%         nszData = openNSx(sources(s).nsfile,1,32,starts(p)-round(0.1*ns.fs),total+round(0.1*ns.fs),'p')';
        nsxData = openNSx(sources(s).nsfile, 'read', 'c:1:32', ['t:', num2str(starts(p)-round(0.1*Fs)), ':', num2str(starts(p)+total+round(0.1*Fs))]);

        tmpDataLen = size(nsxData.Data,2);
        if tmpDataLen ~= total
            nsxData.Data = nsxData.Data(:, 1:(total));
        end
        
        tmpdata = double(nsxData.Data');
%         tmpdata = zeros(total/(Fs/newfs),32,length(allpts));
        % downsample to 5k
        Fs = Fs;  % Sampling Frequency
        
        Fpass = 4000;    % Passband Frequency
        Fstop = 5000;    % Stopband Frequency
        Apass = 1;       % Passband Ripple (dB)
        Astop = 60;      % Stopband Attenuation (dB)
        match = 'both';  % Band to match exactly

        parfor k=1:size(tmpdata,2)
            % Construct an FDESIGN object and call its ELLIP method.
            h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, Fs);
            Hd = design(h, 'ellip', 'MatchExactly', match);
            tmpdata(:,k)=filtfilt(Hd.sosMatrix,Hd.ScaleValues,tmpdata(:,k));
        end
        tmpdata=tmpdata((ns.fs/newFs)/2:(ns.fs/newFs):end,:);

        % define parameters
        bandstop=[... % remove 60 hz noise
            59.2 60.8;...
            119.2 120.8;...
            179.2 180.8;...
            239.2 240.8;...
            299.2 300.8;...
            359.2 360.8;...
            419.2 420.8;...
            479.2 480.8;];
        
        % filtering line noise
        for b=1:size(bandstop,1)
            parfor k=1:size(tmpdata,2)
                
                

                
                Fpass1 = bandstop(b,1)-1; % First Passband Frequency
                Fstop1 = bandstop(b,1);   % First Stopband Frequency
                Fstop2 = bandstop(b,2);   % Second Stopband Frequency
                Fpass2 = bandstop(b,2)+1; % Second Passband Frequency
                Apass1 = 0.5;     % First Passband Ripple (dB)
                Astop  = 60;      % Stopband Attenuation (dB)
                Apass2 = 1;       % Second Passband Ripple (dB)
                match  = 'both';  % Band to match exactly
                
                % Construct an FDESIGN object and call its ELLIP method.
                h  = fdesign.bandstop(Fpass1, Fstop1, Fstop2, Fpass2, Apass1, Astop, ...
                    Apass2, newFs);
                Hd = design(h, 'ellip', 'MatchExactly', match);
                tmpdata(:,k)=filtfilt(Hd.sosMatrix,Hd.ScaleValues,tmpdata(:,k));
                
%                 Hd=bsfilt(bandstop(b,:),newFs);
%                 tmpdata(:,k)=filtfilt(Hd.sosMatrix,Hd.ScaleValues,tmpdata(:,k));
            end
        end
%         tmpdata([1:round+(0.5*newFs)],:)=[];
%         tmpdata([1:round(0.5*newFs), end-round(0.5*newFs)+1:end],:)=[];

        % CAR
        tmpdata(:,1:16)=tmpdata(:,1:16)-repmat(mean(tmpdata(:,1:16),2),1,16); % facemotor
        tmpdata(:,17:32)=tmpdata(:,17:32)-repmat(mean(tmpdata(:,17:32),2),1,16); % wernickes
        
        % save back
        data{s}(:,:,p)=tmpdata;
        clc
        fprintf('Current Day: %d/%d\n', s, 4)
        fprintf('Trial num: %d/%d\n', p, length(allpts))
        
    end
    class{s}=tmpclass;
end
data=cat(3,data{:});
class=cat(2,class{:});

if true
    disp('hello world')
end


% 
save('E:\data\PLI\delta\PLIOutput\Delta_ProcessedTrialData.mat','data','class', 'Header','-v7.3');