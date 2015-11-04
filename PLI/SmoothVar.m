function [ varOut ] = SmoothVar( data)%, varargin )

% if ~isempty(varargin)
%     samp = varargin{1};
% else
%     samp = 5;
% end

varOut = nan(1,length(data));
samp = 5;
for ii = 1:length(data)
    if ii == 1
        varOut(ii) = var(data(ii));     
    elseif ii == 2
        varOut(ii) = var(data(1:3));
    elseif ii == length(data)-1
        varOut(ii) = var(data(end-2:end));
    elseif ii == length(data)
        varOut(ii) = var(data(end));
    else
        varOut(ii) = var(data(ii-1:ii+2));
    end % END IF
    
end % END FOR
end % END FUNCTION

