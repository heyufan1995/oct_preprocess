function mfprintf(fids,varargin)
% fprintf with multiple file
%   one example usage is to print to a file and print to the command line
%   (fid = 1)

if isempty(fids)
    return
end

for i = 1:length(fids)
    fprintf(fids(i),varargin{:});
end