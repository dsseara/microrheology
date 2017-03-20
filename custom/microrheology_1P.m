% microrheology_1P.m
% 
% Calculates the one point microrheology from outputs from
% Maria Kilfoil's feature finding and tracking algorithms
%
% INPUTS   msdparams : structure
%
%
%       msdparams.
%                 timeint     = Time interval between frames in seconds
%                 totalFrames = Total number of frames
%                 rg_max      = Maximum Rg to be included in um (Inf for no max)
%                 rg_min      = Minimum Rg to be included in um (0 for no min)
%                 maxtime     = Maximum time to be output in logarithmically
%                                spaced msd and tau ([] for all time)
%                 tc          = Critical frame number (optional, pass [] for all)
%                 tcFraction  = What fraction of tc to cut data as before and after
%                 numFOV      = Number of fields of view
%                 dedriftBool = true to dedrift, false to not
%                 smoo        = If dedriftBool, number of frames to average drift over
%
% OUTPUTS  msdtau  : Logarithmically spaced msd and lag times
%
% You can find all the software and the instructions that this script follows to a T at:
% http://people.umass.edu/kilfoil/downloads.html
%
% Created by Daniel Seara at 2017/03/05 19:07

%function [msdtau, msd, tau] = microrheology_1P(basepath, msdparams)
cd(basepath);

if isempty(msdparams.maxtime)
    msdparams.maxtime = msdparams.totalFrames * msdparams.timeint;
end
if ~isfield(msdparams,'number_of_frames')
    msdparams.number_of_frames = msdparams.totalFrames;
end

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
tc = floor(msdparams.tcFraction*msdparams.tc);

disp('Calculating MSD and taus')
[msd, msdx, msdy, tau, beadcount] = Mean_SD_many_single_beads(basepath, msdparams.timeint, msdparams.number_of_frames, [msdparams.rg_min, msdparams.rg_max], tc);
disp('Done.');

% if sum(~isnan(MSD))==0
%     error('Got all NaNs for MSD');
% end

% Now space MSD out logarithmically with ~16 points per decade
disp('Spacing MSD logarithmically...')
if isempty(msdparams.tc)
    [msdtau] = making_logarithmically_spaced_msd_vs_tau(msd,tau,msdparams.maxtime);
else
    [msdtau.pre] = making_logarithmically_spaced_msd_vs_tau(msd.pre,tau.pre, msdparams.timeint*(tc-1));
    [msdtau.post] = making_logarithmically_spaced_msd_vs_tau(msd.post,tau.post, (msdparams.maxtime - msdparams.timeint*tc) );
end
disp('Done.')

saveTo = [basepath '1pt_msd_of_'  num2str(beadcount) '_beads' filesep];
[status, message, messageid] = mkdir(saveTo);

if isempty(msdparams.tc)
    save([saveTo 'microrheology_1P.mat'], 'msd','msdx','msdy','tau','msdtau', 'msdparams')
else
    save([saveTo 'tc_' num2str(tc) '_microrheology_1P.mat'], 'msd','msdx','msdy','tau','msdtau', 'msdparams')
end    
