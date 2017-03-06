function[msdtau] = making_logarithmically_spaced_msd_vs_tau(msd,tau,maxtime)
% Converts MSD and tau vectors in one 2 column matrix spaced
% logarithmically in tau, roughly 16 points per decade. Follows 
% "Mean_SD_many_single_beads_p".
%
% INPUTS
% 
% msd - the MSD vector
% tau - the tau vector
% maxtime - the maximum time to consider in second
%
% OUTPUTS
%
% msdtau - a 2 column matrix, with tau in the first column and the MSD in
% the second column.

timestep = min(tau);
dt = round(1.15.^(0:99));
[eldt,idt]=unique(dt);
dt = dt(idt);
w = find( dt <= round(maxtime/timestep) ); ndt=length(w);
if ndt > 0
    dt = dt(w);
else
    warning('Invalid maximum dt!')
end

msdtau(:,1) = tau(dt);
msdtau(:,2) = msd(dt);
