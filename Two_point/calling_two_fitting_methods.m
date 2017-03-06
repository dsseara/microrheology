function[] = calling_two_fitting_methods(basepath, data,tau,minradius,maxradius,displaying,beadradius,savingmsdd,savingplot)
% Computes the 2P MSD from the output of "twopoint" using two methods:
% averaging r*Drr and fitting it to a line (to use the extrapolation to 0
% separation). Also calculates the Poisson ratio from the former method. 
% Should be called once for each tau one is interested in, in a loop; the 
% final result will contain data for all previous taus and the new ones 
% being added one at a time.
% Reinitializes if tau==1.
%
% This is a front end for "fitting_automater_fit_to_line_shell" and
% "fitting_automater_r_times_Drr_Dqq_shell", which are themselves front
% ends for "msddE". The front ends are used to remove points with negative
% correlation from the r*Drr fits. 
% Follows "twopoint"
%
% INPUTS
%
% basepath - the base path for the experiment
% data - output from calling "twopoint"
% tau - index of the lag time tau for which the MSD should be calculated.
% minradius - minimum radius in micrometers which will be used to compute
% the MSD
% maxradius - maximum radius in micrometers which will be used to compute
% the MSD
% displaying - a switch which will allow printing of the MSD data on screen
% if set to 1 (typicallly set to 0)
% beadradius - the probe beads radius, in micrometers.
% savingmsdd - a switch which will allow recording of the MSD in a file if
% set to 1 (typically set to 1)
% savingplot - a switch which will allow saving of r*Drr vs r (and such)
% plots in subfolders if set to 1 (typically set to 1)
%
% OUTPUTS
%
% Creates a series of subfolders starting at "rDrr_rDqq_figs". If "savingplot"
% is set to 1, the subfolders will contain plots of r*Drr vs r,
% r*Dthetatheta vs r for the "averaging r*Drr" method, and r*Drr vs r for the
% "fit to line" method. It will also plot the number of points used vs r
% for each tau. Error bars on r*Drr plots come from the standard deviation
% calculated by "twopoint".
% If savingmsdd is set to 1, it will save a file called "msd2P" which
% contains the msd2P matrix, the flor2P matrix and the pois matrix.
%
% msd2P matrix format:
% msd2P(:,:,1) is the MSD data from the average r*Drr method
% msd2P(:,:,2) is the MSD data from the fit r*Drr to a line method
% Each row corresponds to one tau value
% columns:
% 1 - tau in seconds
% 2:3 - raw MSD from r*Drr and r*Dtheta_theta
% 4:5 - smoothed MSD using lfit on the previous 2 columns (done in this
% file)
% 6   - error estimate from the average standard deviations of the r*Drr
% points
% 7:8 - error estimates for r*Drr and r*Dtheta_theta MSDs based on the
% standard deviation of r*Drr with respect to the average (0 for the fit to
% line method)
%
% pois matrix format:
% one row per tau
% columns:
% 1 - tau in seconds
% 2 - Poisson ratio calculated from the ratio of the raw MSD_rr and
% MSD_theta_theta
% 3 - Poisson ratio calculated from the ratio of the smoothed MSD_rr and
% MSD_theta_theta
% 
% flor2P matrix format: (directly output from msddE)
% one row per tau
% columns:
% 1:2 - from fit to line method, slope*(2/(3*beadradius))*(average r) for
% Drr and Dtheta_theta

width=0.9; % for smoothing the msd


path=[basepath 'rDrr_rDqq_figs\'];

rDrrpath = [path 'rDrr_plots\'];
linfitpath = [path 'rDrr_fit_to_line'];

parttitle = ['max radius = ' num2str(maxradius) '\mum'];

[status, message, messageid] = mkdir(rDrrpath);
[status, message, messageid] = mkdir(linfitpath);
% [status, message, messageid] = mkdir([linfitpath '\lin_lin_plots']);
[status, message, messageid] = mkdir([path '\Number_of_points_vs_radius']);



[msdd1 flor1] = fitting_automater_fit_to_line_shell(data,tau,minradius,maxradius,displaying,beadradius,linfitpath,savingplot);

[msdd0 flor0] = fitting_automater_r_times_Drr_Dqq_shell(data,tau,minradius,maxradius,displaying,beadradius,rDrrpath,savingplot);

number_of_points_per_r(data,tau,path,parttitle)

if savingmsdd == 1
    if tau ~= 1
        load([path 'msd2P'])
    end
    msd2P(tau,:,1)=msdd0;       % Average of r*Drr
    msd2P(tau,:,2)=msdd1;       % Fit r*Drr to a line
warning( 'off', 'MATLAB:polyfit:PolyNotUnique' );
    msd2P(:,4,1) = lfit(msd2P(:,1,1),msd2P(:,2,1),width);
    msd2P(:,5,1) = lfit(msd2P(:,1,1),msd2P(:,3,1),width);
%     msd2P(:,4,2) = lfit(msd2P(:,1,2),msd2P(:,2,2),width);     % commented out, too many  warning messages!
%     msd2P(:,5,2) = lfit(msd2P(:,1,2),msd2P(:,3,2),width);
warning( 'on', 'MATLAB:polyfit:PolyNotUnique' );
    flor2P(tau,:,1)=flor0;
    flor2P(tau,:,2)=flor1;
    
    pois=poisson_from_msdd( msd2P );
    save([path 'msd2P.mat'],'msd2P','flor2P', 'pois')
end