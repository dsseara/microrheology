% microrheology_1P.m
% 
% Calculates the one point microrheology from outputs from
% Maria Kilfoil's feature finding and tracking algorithms
%
% INPUTS   msdparams : struct with parameters to be used in 1 point MSD calculation
%
%          FIELD
%          timeint          = Time interval between frames in seconds
%          number_of_frames = Total number of frames
%          rg_cutoff        = Maximum Rg to be included in um ([] for all)
%          maxtime          = Maximum time to be output (in logarithmically spaced msd and tau
%          tc               = Critical frame number (optional, pass [] for entire series)
%          numFOV           = Number of fields of view
%          dedriftBool      = true to dedrift, false to not
%          smoo             = If dedriftBool, number of frames to average drift over 
%
% OUTPUTS  msdtau  : Logarithmically spaced msd and lag times
%
% You can find all the software and the instructions that this script follows to a T at:
% http://people.umass.edu/kilfoil/downloads.html
%
% Created by Daniel Seara at 2017/03/05 19:07

%function [msdtau, msd, tau] = microrheology_1P(basepath, msdparams)
cd(basepath);

converted = 0;
while converted ==0
    converted = input('Need to hard code pixels-to-micron conversion in pixtomicro.m. \n When done, enter 1: ');
end

disp('Dedrifting and converting...')
if msdparams.dedriftBool
    for jj = 1:msdparams.numFOV
        dedrifting_and_conversions(basepath,jj,smoo);
    end
else
    for jj=1:msdparams.numFOV
        conversions_no_dd(basepath,jj);
    end
end
disp('Done.')

%%
% Separate data by individual bead tracks

disp('Sorting by tracks...')
getting_individual_beads(basepath,1:msdparams.numFOV);
disp('Done')

%%
% Now calculate the radius of gyration for each bead

disp('Calculating radius of gyration for beads...')
rg_matrix_many_single_beads(basepath);
disp('Done.');

%%
% Now actually get the 1P MSD, computed as sum of x and y MSDs. Will then be converted to complex viscoelastic modulus

disp('Calculating MSD and taus')
[msd, tau] = Mean_SD_many_single_beads( basepath, msdparams.timeint, msdparams.number_of_frames, msdparams.rg_cutoff, floor((3/3)*msdparams.tc) );
disp('Done.');


if ispc
    save([basepath, '\1pt_msd\parameters.mat'], 'msdparams')
elseif isunix
    save([basepath, '/1pt_msd/parameters.mat'], 'msdparams')
end

% if sum(~isnan(MSD))==0
%     error('Got all NaNs for MSD');
% end

%%
% Now space MSD out logarithmically with ~16 points per decade

disp('Spacing MSD logarithmically...')
if isempty(msdparams.tc)
    [msdtau] = making_logarithmically_spaced_msd_vs_tau(msd,tau,msdparams.maxtime);
else
    [msdtau.pre] = making_logarithmically_spaced_msd_vs_tau(msd.pre,tau.pre,5*(floor((3/3)*msdparams.tc)));
    [msdtau.post] = making_logarithmically_spaced_msd_vs_tau(msd.post,tau.post,(msdparams.maxtime - 5*(floor((3/3)*msdparams.tc))));
end
disp('Done.')

if ispc
    save([basepath, '\1pt_msd\logarithmically_spaced_msd.mat'], 'msdtau')
elseif isunix
    save([basepath, '/1pt_msd/logarithmically_spaced_msd.mat'], 'msdtau')
end
%end

