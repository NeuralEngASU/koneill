%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cohereEpilepsy
%   Desc: Calculates the magnitude squared coherence between 
%
%   [P,R] = PLI(D1,D2)
%   For two time-domain signals D1 and D2, calculate the phase-lag index P
%   and the phase-similarity measure R.  D1 and D2 must have the same
%   number of columns; P and R will be vectors with the same number of
%   entries as the number of columns in D1 and D2.
%
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [c] = cohereEpilepsy(data1,data2, cParams)

        Cxy =  mscohere(data1, data2,hamming(cParams.sizeWindow),cParams.numOverlap,cParams.freqBand,cParams.Fs);
        c = mean(Cxy)';
    
%     chronuxParams.Fs = 500;
%     chronuxParams.fpass = [cParams.freqBand(1), cParams.freqBand(end)];
%     chronuxParams.pad = 0;
%     chronuxParams.tapers = [3, 5];
%     [~,Cmat,~,~,~,f]=CrossSpecMatc([data1, data2],1,chronuxParams);
%         
%     fIdx = f>=cParams.freqBand(1) & f<=cParams.freqBand(end);
%     
%     c = mean(abs(Cmat(fIdx,1,2)))';
end % END FUNCTION cohereEpilepsy

% EOF