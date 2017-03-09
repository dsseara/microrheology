function pois = poisson_from_msdd(msd2P)
%
%       A little function for calculating the Poisson ratio out of MSD_par,MSD_perp 
%

% a useful ratio
ratio=msd2P(:,3,1)./2./msd2P(:,2,1);
ratio_smooth=msd2P(:,5,1)./2./msd2P(:,4,1);

pois(:,1)=msd2P(:,1,1);
pois(:,2)=(4*ratio-3)./(4*ratio-4);
pois(:,3)=(4*ratio_smooth-3)./(4*ratio_smooth-4);

end
