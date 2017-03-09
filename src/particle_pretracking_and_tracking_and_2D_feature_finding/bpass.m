function s = bpass(img,lambda,w)
%      7-23-03  Maria Kilfoil
%
% 		Implements a real-space bandpass filter to suppress pixel noise and 
%       slow-scale image variations while retaining information of a characteristic size.
% 		*Works with anisotropic 3d cube data*
% 
%  CALLING SEQUENCE:
% 		res = bpass(img, lambda, w)
%  INPUTS:
% 		img:	two-dimensional array to be filtered.
% 		lambda: characteristic lengthscale of noise in pixels. Additive noise averaged 
%                   over this length should vanish. May assume any positive floating value.
% 			Make it a 3-vector if aspect ratio is not 1:1:1.
% 		w: A length in pixels somewhat larger than *half* a typical object. Must be an odd valued 
%                   integer. Make it a 3-vector if aspect ratio is not 1:1:1.
%  OUTPUTS:
% 		res:	filtered image.
%  PROCEDURE:
% 		simple 'mexican hat' wavelet convolution yields spatial bandpass filtering.
%  NOTES:
%		based on "Methods of digital video microscopy for colloidal studies", John Crocker 
%       and David Grier, J. Colloid Interface Sci. 179, 298 (1996), and on bpass.pro IDL code 
%       written by John Crocker and David Grier. 
%
clear s t
a=double(img);
b=double(lambda);
w = round(max(w,2*b));
N = 2*w + 1;
r = [-w:w]/(2*b);
xpt = exp(-r.^2);  
B=(sum(xpt))^2;
xpt = xpt / sum(xpt);    %xpt=exp(-[-w:w].^2/4*lambda^2)/sqrt(B)
factor=((sum(xpt.^2))^2-1/N^2)*B;  %sum(xpt.^2)=1/Bexp(-(i^2+j^2)/(4*lambda^2))
% note: N not N^2 etc since doing 2D conv along each axis separately
gx = xpt;    
gy=gx';
bx = zeros(1,N)-1./N;  
by=bx';
% g = conv2( a, gx );
% g = conv2( g, gy );
% b = conv2( a, bx );
% b = conv2( b, by );
g = conv2(gx,gy,a,'valid');
b = conv2(bx,by,a,'valid');
res = g-b;
s=max(res/factor,0);
% whos s
tmp=zeros(length(a(:,1)),length(a(1,:)));
tmp(w+1:(length(s(:,1))+w),w+1:(length(s(1,:))+w))=s;
s=tmp; %filtered image with w width zero border
