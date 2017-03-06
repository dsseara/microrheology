% This script calculates the one point microrheology from outputs from Maria Kilfoil's feature finding and tracking algorithms
%
% You can find all the software and the instructions that this script follows to a T at:
% http://people.umass.edu/kilfoil/downloads.html
%
% Created by Daniel Seara at 2017/03/05 19:07

basepath = input('Enter the path to the experimental data \n (make sure it ends with the appropriate slash, pwd won"t work): ');
cd(basepath);

dedrift_bool = input('Enter 1 or 0 if data does or does not need to be dedrifted: ');
numFOV = input('numFOV: Enter the number of fields of view for the data: ');

converted = 0;

while converted ==0
    converted = input('Need to hard code pixels-to-micron conversion in pixtomicro.m. \n When done, enter 1: ');
end

disp('Dedrifting and converting...')
if dedrift_bool
    smoo = input('smoo: Number of frames the drift is averaged before subtracted (should be pretty large): ');
    for jj = 1:numFOV
        dedrifting_and_conversions(basepath,jj,smoo);
    end
else
    for jj=1:numFOV
        conversions_no_dd(basepath,jj);
    end
end
disp('Done.')

%%
% Separate data by individual bead tracks

disp('Sorting by tracks...')
getting_individual_beads(basepath,1:numFOV);
disp('Done')

%%
% Now calculate the radius of gyration for each bead

disp('Calculating radius of gyration for beads...')
rg_matrix_many_single_beads(basepath);
disp('Done.');

%%
% Now actually get the 1P MSD, computed as sum of x and y MSDs. Will then be converted to complex viscoelastic modulus
disp('Enter 1P MSD parameters...')
msdparams.timeint          = input('timeint: Time interval between frames in seconds: ');
msdparams.number_of_frames = input('number_of_frames: Maximum lag time (typically total number of frames available): ');
msdparams.rg_cutoff        = input('rg_cutoff: Maximum Rg to be included in um (leave empty or set to negative number if unneeded): ');
msdparams.maxtime          = input('maxtime: Maximum time to be output (in logarithmically spaced msd and tau): ');

disp('Calculating MSD and taus')
[MSD, tau] = Mean_SD_many_single_beads( basepath, msdparams.timeint, msdparams.number_of_frames, msdparams.rg_cutoff );
disp('Done.');

if ispc
    save([basepath, '\1pt_msd\parameters.mat'], 'msdparams')
elseif isunix
    save([basepath, '/1pt_msd/parameters.mat'], 'msdparams')
end

if sum(~isnan(MSD))==0
    error('Got all NaNs for MSD');
end

%%
% No space MSD out logarithmically with ~16 points per decade

disp('Spacing MSD logarithmically...')
[msdtau] = making_logarithmically_spaced_msd_vs_tau(MSD,tau,msdparams.maxtime);
disp('Done.')



