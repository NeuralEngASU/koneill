sources(1).marker='E:\data\human CNS\Utah\delta\verbal\MARKERS_day1.mat';
sources(1).nsfile='E:\data\human CNS\Utah\delta\verbal\20080730-151259\20080730-151259-002.ns5';
sources(2).marker='E:\data\human CNS\Utah\delta\verbal\MARKERS_day2.mat';
sources(2).nsfile='E:\data\human CNS\Utah\delta\verbal\20080731-111811\20080731-111811-001.ns5';
sources(3).marker='E:\data\human CNS\Utah\delta\verbal\MARKERS_day3.mat';
sources(3).nsfile='E:\data\human CNS\Utah\delta\verbal\20080801-102946\20080801-102946-001.ns5';
sources(4).marker='E:\data\human CNS\Utah\delta\verbal\MARKERS_day4.mat';
sources(4).nsfile='E:\data\human CNS\Utah\delta\verbal\20080801-135921\20080801-135921-001.ns5';

for s = 1:4
    
    load(sources(s).marker,'offset_words','pts_words','lbl_words');
    data = openNSx(sources(s).nsfile, 'read', 'c:1:32');
    Header = data.MetaTags;
    data = data.Data;
    
    save(['D:\PLI\Speech\DeltaSpeech_Day', num2str(s), '.mat'], 'data', 'Header', 'lbl_words', 'offset_words', 'pts_words', '-v7.3')
    
end % END FOR



% EOF