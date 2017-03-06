function r = localmax_nodilate(image,sep,field,intmin)
%	identical to John's version of local_max
range = double(fix(sep/2));
% whos intmin
% intmin
clear a h
% a=int16((im-min(min(im)))*255/(max(max(im))-min(min(im))));
a=image;
w = round( 2 * range + 1 );
s = rsqd( w,w );
mask=(s <= range^2);
yrange = range;
if (or(field==1,field==0))==1
    mask = fieldof( mask, 0 );
    yrange = double(int16(range/2.)) +1;
elseif field == 2
    ;
else
    error('Field parameter in localmax must be 0, 1 or 2');
end
% b = imdilate( a, mask);
clear q
% but don't include pixels from the
% background which will be too dim
if intmin== 0
%     whos a
% else
%   ;
    h=sum(hist(double(a),256),2);  %a is the filtered image
    h=cumsum(h);
    h=h/max(h);
    intmin = 1 + sum(h < 0.70);
%     while (h(intmin) < 0.64),
%         intmin = intmin + 1;
%     end
%   ** this may not be robust; forces min intensity to be some fraction of the total
%   (more precisely, throws out bottom 1/5 of pixels
    if (intmin < 5) 
        intmin=round(length(h)/5);
    end
end
%intmin
% r = find(a==b & a >= intmin);
r = find(a >= intmin);
% size(r)

% Discard maxima within range of the edge
sz = size( a );
nx = sz(2);
ny = sz(1);
y = mod(r,ny);
x = fix(r / ny+1);
x0 = double(x) - range;
x1 = double(x) + range;
y0 = double(y) - yrange;
y1 = double(y) + yrange;
good = find( (x0 >= 1) & (x1 < nx) & (y0 >= 1) & (y1 < ny));
ngood=length(good);
r = r(good);
x = x(good);
y = y(good);
x0 = x0(good);
x1 = x1(good);
y0 = y0(good);
y1 = y1(good);
% Find and clear spurious points arising from features which get 
% found twice or which have flat peaks and thus produce multiple hits.  
c=zeros(ny,nx);
c(r) = a(r);
center = w * range + range + 1; % position in mask pixel number of the center of the mask
for i = 1:length(r),
    b = c(y0(i):y1(i),x0(i):x1(i));
    b =  b.*mask;  % 	look only in circular region
    [Y,I] = sort(b,1);
    g1=I(length(I(:,1)),:); % array containing in which row was the biggest value for each column
    [yi,locx]=max(Y(length(Y(:,1)),:));  % yi is abs maximum within mask, locx is column #
    % *** come back to fix this, since will lose local maxima that are in the
    % center, if they are not larger than a fluctuation nearby (within the
    % radius) ***
    locy = g1(locx);  % row number corresponding to maximum
    [d1 d2]=size(b);
    location = (locx-1)*d1+locy;  % find location of remaining maximum in the mask
    if (location ~= center), %compare location of max to position of center of max
		c(y(i),x(i)) = 0;
    end
end
r = find( c ~= 0 );	% What's left are valid maxima.
% return their locations
% find coords x = uint16(r / ny+1);
% y = mod(r,ny);
% size(r)