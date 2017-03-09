function theta = thetarr(w)
% produce a theta mask
theta = zeros(w,w);
xc = (w-1)/2;
yc = (w-1)/2;
x = [-xc:xc];
y = [-yc:yc];
% for j = 0, w-1 do begin  ** switched array start to 1
for j = 1:w,
	theta(:,j) = double(atan2(x',y(j)));
end
