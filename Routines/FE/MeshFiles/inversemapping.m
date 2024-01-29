function map = inversemapping(idx,src,N)
if exist('N','var')==0 || isempty(idx)
    N = length(idx);
end
map = nan(N,1);
map(idx) = 1:nnz(idx);

if exist('src','var') && ~isempty(src)
    map = map(src);
end
end