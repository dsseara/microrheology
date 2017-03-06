function r = feature2D(img,lambda,w,masscut,Imin,field)
% Note: feature2D requires the Image processing toolbox due to a call to
% imdilate in the localmax subfunction. Use feature2D_nodilate for an
% alternative which works almost as well.

%      7-29-03  Maria Kilfoil
%extent should be 2*w+1 in which w is the same as bpass.
% 	Finds and measures roughly circular 'features' within an image.
%  CALLING SEQUENCE:
% 	f = feature2D(image,lambda,diameter,masscut,Imin,field)
%  INPUTS:
% 	img:	(nx,ny) array which presumably contains some features worth finding
%   lambda: length scale of noise to be filtered out, in pixels; typically 1
% 	w:      a parameter which should be a little greater than the radius of the 
%           largest features in the image.
% 	Imin: 	(optional) Set this optional parameter to the minimum allowed value for the peak 
%           brightness of a feature. Useful for limiting the number of spurious features in
%       	noisy images. If set to 0 (default), the top 30% of the bright pixels
%       	will be used.
% 	masscut: (optional) Setting this parameter saves runtime by reducing the runtime wasted on 
%           low mass 'noise' features. (default is 0, accept all.)
% 	field: 	(optional) Set this parameter to 0 or 1 if image is actually just one (odd or even) field of an interlaced 
%           (e.g. video) image. All the masks will then be constructed with a 2:1 aspect ratio. Otherwise 
%           set to 2 for progressive scan cameras. 
%
% NOT IMPLEMENTED:
% 	separation: an optional parameter which specifies the minimum allowable separation 
%           between feature centers. The default value is diameter+1.
%
%
%  OUTPUTS:
% 		f(:,1):	the x centroid positions, in pixels.
% 		f(:,2): the y centroid positions, in pixels. 
% 		f(:,3): integrated brightness of the features. ("mass")
% 		f(:,4): the square of the radius of gyration of the features.
% 		    (second moment of the "mass" distribution, where mass=intensity)
% 		f(:,5): eccentricity, which should be zero for circularly symmetric features and 
%                   order one for very elongated images.
%  RESTRICTIONS:
%       To work properly, the image must consist of bright, circularly symmetric regions 
%       on a roughly zero-valued background. To find dark features, the image should be 
%       inverted and the background subtracted. If the image contains a large amount of 
%       high spatial frequency noise, performance will be improved by first filtering the image.
%       BPASS will remove high spatial frequency noise, and subtract the image background. 
%       Individual features should NOT overlap.
%
%  MODIFICATION HISTORY:
% 		This code is inspired by feature_stats2 written by
% 			David G. Grier, U of Chicago, 			 1992.
% 		Written by John C. Crocker, U of Chicago, optimizing 
% 			runtime and measurement error, 			10/93.
%       	Matlab version written by Maria L. Kilfoil		2003.
%       10-18-07 Vincent Pelletier, Maria Kilfoil -- added masscut to Matlab version

if nargin < 6, field=2;, end
if nargin < 5, Imin=0;, end
if nargin < 4, masscut=0;, end

extent=2*w+1;
image = bpass(img,lambda,w);
if (mod(extent,2) == 0),
    disp('Requires an odd extent.  Adding 1...');
    extent = extent + 1;
end
sz = size(image);
nx = sz(2);
ny = sz(1);
% if n_params() eq 2 then sep = extent+1
sep = extent; 

%       Put a border around the image to prevent mask out-of-bounds
% Only use the following 2 lines if you are not using bpass to do spatial filtering first.
% Otherwise, image returned from bpass already had a border of w width
% a = zeros( ny + extent, nx + extent );
% a(fix(extent/2)+1:(fix(extent/2))+ny,fix(extent/2)+1:(fix(extent/2))+nx) = image;
a=image;

%       Finding local maxima
loc= localmax(image,sep,field,Imin);
if (numel(loc) == 0 || loc(1) == -1)
    r = -1;
    return
else
    ;
end
y = mod(loc,ny);
x = fix(loc / ny+1);  %the pixel position of the maximum
nmax=length(loc);         
m=zeros(length(loc),1);
xl = x - fix(extent/2);
xh = xl + extent -1;        

%       Set up some masks
rsq = rsqd( extent,extent );
t = thetarr( extent );

mask = le(rsq,(extent/2)^2);      
mask2 = ones(1,extent)'*[1:extent];
mask2 = mask2.*mask;           
mask3= (rsq.*mask) + (1/6);
cen = (extent-1)/2 +1;           
% cmask = vpa(cos(sym('2')*t)).*mask;  
% ultra high presision, since Matlab and IDL differ here
cmask = cos(2*t).*mask;
% smask = vpa(sin(sym('2')*t)).*mask;
smask = sin(2*t).*mask;
cmask(cen,cen) = 0.0;
smask(cen,cen) = 0.0;  

suba = zeros(extent, extent, nmax);
xmask = mask2;
ymask = mask2';
yl = y - fix(extent/2);
yh = yl + extent -1;           
yscale = 1;
ycen = cen;                   

%	Estimate the mass	
for i=1:nmax, 
    m(i) = sum(sum(double(a(fix(yl(i)):fix(yh(i)),fix(xl(i)):fix(xh(i)))).*mask)); 
end

% remove features based on 'masscut' parameter
b = find(m > masscut);    %only those features with a total mass higher than masscut make the grade
nmax=length(b);
if nmax==0 
    disp('No feature found!');
    r=[];
    return
end
xl = xl(b);
xh = xh(b);
yl = yl(b);
yh = yh(b);
x = x(b);
y = y(b);
m = m(b);

disp(strcat(num2str(nmax,'%01.0f'),' features found.'));

%	Setup some result arrays
xc = zeros(nmax,1);
yc = zeros(nmax,1);
rg = zeros(nmax,1);
e  = zeros(nmax,1);
%	Calculate feature centers
for i=1:nmax,
	xc(i) = sum(sum(double(a(fix(yl(i)):fix(yh(i)),fix(xl(i)):fix(xh(i)))).*xmask));  
	yc(i) = sum(sum(double(a(fix(yl(i)):fix(yh(i)),fix(xl(i)):fix(xh(i)))).*ymask));
end
x1=x;
y1=y;
%	Correct for the 'offset' of the centroid masks
xc = xc./m - ((extent+1)/2);             
yc = (yc./m - (extent+1)/2)/yscale;  
%	Update the positions and correct for the width of the 'border'
x = x + xc - 0*fix(extent/2);
y = ( y + yc - 0*fix(extent/2) ) * yscale;
x2=x;
y2=y;

%	Construct the subarray and calculate the mass, squared radius of gyration, eccentricity
for i=1:nmax,
    suba(:,:,i) = fracshift( double(a(fix(yl(i)):fix(yh(i)),fix(xl(i)):fix(xh(i)))), -xc(i) , -yc(i) );
    m(i) = sum(sum(( suba(:,:,i).*mask )));             % mass
    rg(i) = (sum(sum( suba(:,:,i).*mask3 ))) / m(i);    % squared radius of gyration
    tmp = sqrt(( (sum(sum( suba(:,:,i).*cmask )))^2 ) +( (sum(sum( suba(:,:,i).*smask )))^2 )); 
    tmp2 = (m(i)-suba(cen,ycen,i)+1e-6);
    e(i) = tmp/tmp2;                                    % eccentricity
end
for i=1:nmax,
	xc(i) = sum(sum(double(suba(:,:,i)).*xmask));  
	yc(i) = sum(sum(double(suba(:,:,i)).*ymask));
end
xc = xc./m - ((extent+1)/2);             
yc = (yc./m - (extent+1)/2)/yscale;  %get mass center
x3 = x2 + xc - 0*fix(extent/2);
y3 = ( y2 + yc - 0*fix(extent/2) ) * yscale;

r = [x3,y3,m,rg,e];
%toc