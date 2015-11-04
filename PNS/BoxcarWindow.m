function xx = BoxcarWindow( len, fs )
% Causal Boxcar Window
%   length (ms) , fs
%
%              _________
% ____________|         |________
%             0        length
%

% xx = zeros(1, fs/1000*len+2);
xx = zeros(1, fs/1000*len+1);

% xx(1) = 0;
% xx(2:end-1) = 1;
% xx(end) = 0;

xx(1) = 0;
xx(2:end) = 1;
% xx(end) = 0;


end

