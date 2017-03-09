function[resmmm] = pixtomicro(res)

%For our purposes, the plots should be in micrometers and not in pixels.
%The following give the conversion factors (as written on microscope used
%in Summer 2006):
%10X ===> 1 pixel = (0.954 +/- 0.048) micrometers
%40X ===> 1 pixel = (0.260 +/- 0.018) micrometers
%63X ===> 1 pixel = (0.158 +/- 0.05) micrometers
%This program makes the appropriate conversion. Note that it
%will not take the precision into account. This should be done otherwise.
%Note that this assumes that we took the pictures with a full screen
%resolution, which may or may not be the case. Also if there is 2x2
%binning, as is often the case, there needs to be an additional factor of 2
%in front.
%
% Murrell Lab conversion: 1 pixel = 0.108 microns
resmmm = res;
resmmm(:,1:2) = 2*0.108*res(:,1:2);

end