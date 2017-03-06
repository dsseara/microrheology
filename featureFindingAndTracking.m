% This script runs Maria Kilfoil's Matlab tools for feature finding, particle tracking
%
% You can find all the software and the instructions that this script follows to a T at:
% http://people.umass.edu/kilfoil/downloads.html
%
% Created by Daniel Seara at 2017/03/05 16:59

basepath = input('Enter the path to the experimental data \n (make sure it ends with the appropriate slash, pwd won"t work): ');
cd(basepath);
fprintf('\n Within the path you entered above, there should be folders called "fov#", \n and the frames of the movie are expected to be named "fov#_####.tif". \n' )
numFOV = input('Enter the number of fields of view for the data: ');
for i=1:numFOV
    time = input('Enter a vector of times for the frames: ');
    save(sprintf('fov%i_times',i),'time');
end

%%
% Now input the parameters to be used for feature finding \n in a struct called ffparams (we will be editting these until we get something acceptable). The first seven are editted by hand before moving on.

disp('Enter feature finding parameters...')

ffparams.featsize  = input('featsize: Radius of features in pixels (integer): ');
ffparams.barint    = input('barint: Minimum integrated intensity to be accepted: ');
ffparams.barrg     = input('barrg: Maximum radius of gyration squared to be accepted in pixels squared: ');
ffparams.barcc     = input('barcc: Maximum eccentricity to be accepted: ');
ffparams.IdivRg    = input('IdivRg: Minimum ratio of integrated intensity to radius of gyration squared to be accepted: ');
ffparams.Imin      = input('Imin: Minimum intensity of local maximum to be considered (set to 0 for default of top 30%): ');
ffparams.masscut   = input('masscut: Parameter which defines threshold for integrated intensity before position refinement (optional, set to 0 to skip): ');

ffparams.fovn      = input('fovn: Which field of view are you using to optimize feature finding: ');
ffparams.numframes = input('numframes: Total number of frames: ');
ffparams.frame     = input('frame: Which frame in that fov you are using to optimize (i.e., if fovn = 1 and frame = 3, will look at image fov1/fov1_0003.tif): ');
ffparams.field     = input('field: Enter 0 or 1 if the frame is one field of an interlaced video frame, and 2 if it is a full frame: ');

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
for ii=1:numFOV
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

disp('Enter tracking parameters')
trackparams.featsize   = ffparams.featsize;
trackparams.maxdisp    = input('maxdisp: Maximum displacement in pixels a feature may undergo between successive frames (default 2): ');
trackparams.goodenough = input('goodenough: Minimum length for a trajectory to be retained (default 100): ');
trackparams.memory     = input('memory: Maximum number of frames a feature is allowed to skip (default 1): ');

disp('Tracking...')
for jj = 1:numFOV
    fancytrack(basepath, jj, trackparams.featsize, trackparams.maxdisp, trackparams.goodenough, trackparams.memory );
end
disp('Done.')

if ispc
    save([basepath, 'Bead_tracking\parameters.mat'],'trackparams')
elseif isunix
    save([basepath, 'Bead_tracking/parameters.mat'],'trackparams')
end

fprintf('\n Make sure that the resulting files aren"t empty... \n')