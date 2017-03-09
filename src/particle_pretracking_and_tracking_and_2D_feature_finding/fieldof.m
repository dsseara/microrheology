function res = fieldof(array,f)
% returns the even or odd field of an image
% f=0 for odd fields starting with first, f=1 for even fields
nd = ndims(array);
sz = size(array);
if nd ~= 2
    error('Argument must be a two-dimensional array!')
end
% if keyword_set(odd) then f=1 else f=0
ny2 = uint16( (sz(1)+(1-f))/2 );
rows = [0:double(ny2)-1]*2 +1 + f;
res=array(rows,:);
