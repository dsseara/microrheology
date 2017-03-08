% This script runs Maria Kilfoil's Matlab tools for feature finding, particle tracking
%
% Requires two structs for feature finding (ffparams) and for tracking (trackparams)
% 
%     ffparams. 
%              featsize  = Radius of features in pixels (integer)
%              barint    = Minimum integrated intensity to be accepted
%              barrg     = Maximum radius of gyration squared to be accepted in pixels squared
%              barcc     = Maximum eccentricity to be accepted
%              IdivRg    = Minimum ratio of integrated intensity to radius of gyration squared to be accepted
%              Imin      = Minimum intensity of local maximum to be considered (set to 0 for default of top 30%)
%              masscut   = Parameter which defines threshold for integrated intensity before position refinement (optional, set to skip
%              fovn      = Field of view are you using to optimize feature finding
%              numframes = Total number of frames
%              numFOV    = Number of fields of view for the data
%              frame     = Frame in that fov you are using to optimize (i.e., if fovn = 1 and frame = 3, will look at image fov1_0003.tif)
%              field     = Enter 0 or 1 if the frame is one field of an interlaced video frame, and 2 if it is a full frame
%
% trackparams.
%             featsize   = Same as ffparams.featsize
%             maxdisp    = Maximum displacement in pixels a feature may undergo between successive frames (default 2)
%             goodenough = Minimum length for a trajectory to be retained (default 100)
%             memory     = Maximum number of frames a feature is allowed to skip (default 1)
% You can find all the software and the instructions that this script follows to a T at:
% http://people.umass.edu/kilfoil/downloads.html
%
% Created by Daniel Seara at 2017/03/05 16:59

basepath = input('Enter the path to the experimental data \n (make sure it ends with the appropriate slash, pwd won"t work): ');
cd(basepath);
fprintf('\n Within the path you entered above, there should be folders called "fov#", \n and the frames of the movie are expected to be named "fov#_####.tif". \n' )

for i=1:ffparams.numFOV
    time = input('Enter a vector of times for the frames: ');
    save(sprintf('fov%i_times',i),'time');
end

%%
% Run the feature finding initial step, optimize parameters until satisfied
disp('Open up the image you want to optimize with in ImageJ')
pause(1);

moveOn = 0;
while moveOn == 0
    close all;
    [M2, MT] = mpretrack_init(basepath, ffparams.featsize, ffparams.barint, ffparams.barrg, ffparams.barcc, ffparams.IdivRg, ffparams.fovn, ffparams.frame, ffparams.Imin, ffparams.masscut, ffparams.field);
    moveOn = input('Does the image look good? Enter 0 for no, 1 for yes: ');
    if moveOn==0
        fprintf('\n Change whatever parameters you want, enter "dbcont" when finished \n')
        keyboard
    end
end

disp('Finding features...')
for ii=1:ffparams.numFOV
    mpretrack( basepath, ii, ffparams.featsize, ffparams.barint, ffparams.barrg, ffparams.barcc, ffparams.IdivRg, ffparams.numframes,ffparams.Imin, ffparams.masscut, ffparams.field );
end
disp('Done.')

if ispc
    save([basepath, 'Feature_finding\parameters.mat'],'ffparams')
elseif isunix
    save([basepath, 'Feature_finding/parameters.mat'],'ffparams')
end

%% 
% Now we do the particle tracking, saving the parameters used in a struct called trackparams
disp('Tracking...')
for jj = 1:ffparams.numFOV
    fancytrack(basepath, jj, trackparams.featsize, trackparams.maxdisp, trackparams.goodenough, trackparams.memory );
end
disp('Done.')

if ispc
    save([basepath, 'Bead_tracking\parameters.mat'],'trackparams')
elseif isunix
    save([basepath, 'Bead_tracking/parameters.mat'],'trackparams')
end

fprintf('\n Make sure that the resulting files aren"t empty... \n')