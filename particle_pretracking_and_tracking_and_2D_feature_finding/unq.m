%THIS CODE IS PLANNED FOR ARRAY BEING A ROW VECTOR

function [ret] = unq(array,idx)

s = size(array);
if s, else, warning('array must be an array'), end    %warning if s is empty
if idx
    q = array(idx);
    qshift = circshift(q,[0,-1]);
    indices = find(q~=qshift);
    if indices, ret = idx(indices);, else, ret = length(q);, end
else
    array=array;
    arrayshift = circshift(array,[0,-1]);
    indices = find(array~=arrayshift);
    if indices, ret = indices;, else, ret = length(array);, end
end