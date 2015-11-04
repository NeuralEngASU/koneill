function [p,r,deltaPhi] = wpli(data1,data2)
% PLI Calculate phase-lag index and phase similarity measure
%
%   [P,R] = PLI(D1,D2)
%   For two time-domain signals D1 and D2, calculate the phase-lag index P
%   and the phase-similarity measure R.  D1 and D2 must have the same
%   number of columns; P and R will be vectors with the same number of
%   entries as the number of columns in D1 and D2.

% instantaneous phase
phi1 = atan(imag(hilbert(data1))./data1);
phi2 = atan(imag(hilbert(data2))./data2);

% instantaneous phase difference
deltaPhi = phi2 - phi1;

% phase similarity measure - average length of the phase-difference vector
r = abs(mean(exp(1i*deltaPhi),2));

% phase-lag index - the average sign of the phase difference
% p = abs(mean(sign(deltaPhi),2));

p = abs(mean(abs(deltaPhi).*sign(deltaPhi),2))./mean(abs(deltaPhi),2);