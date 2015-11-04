function data = mat2hdf5( data, varargin )

path = '';

% parse varargin
for i = 1:2:length(varargin)
    eval([lower(varargin{i}) '=varargin{i+1};']);
end

data.streams.(name).data;

dataLen = data.streams.(1).data;

for i = length(data.streams)
    h5create(path, ['/C', num2str(i)], [1, dataLen]);
    h5write(path,  ['/C', num2str(i)], data.streams.(i).data);
end % END FOR

% for i = N:K
%     voltData = h5read(data.info.hdf5Path, ['/C', num2str(i)], timeIdx(1), length(timeIdx));
% end % END FOR

data.info.hdf5Path = fullfile(path, fileName);
data = rmfield(data, ' streams');

end % END FUNCTION

% EOF