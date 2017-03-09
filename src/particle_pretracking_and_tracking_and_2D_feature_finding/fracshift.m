function res = fracshift(im,shiftx,shifty)
%	barrel "shifts" a floating point arr by a fractional pixel amount,
%		by using a 'lego' interpolation technique.
ipx = double(fix( shiftx ));  % integer part
ipy = double(fix( shifty ));
fpx = shiftx - ipx;   % decimal part
fpy = shifty - double(ipy);
% to handle negative shifts:
if fpx < 0
	fpx=fpx+1;
    ipx=ipx-1;
end		
if fpy < 0
	fpy=fpy+1;
    ipy=ipy-1;
end

image=double(im);
imagex  = circshift( image,[ipy ipx+1]   );
imagey  = circshift( image,[ipy+1 ipx]   );
imagexy = circshift( image,[ipy+1 ipx+1] );
image   = circshift( image,[ipy ipx]    );

res = ((1-fpx)*(1-fpy)*image)+ (fpx*(1-fpy)*imagex)+ ((1-fpx)*fpy*imagey)+ (fpx*fpy*imagexy);
