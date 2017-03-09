% This script calculates the Non-gaussian parameter for all taus available
%
% Assumes that featureFindingAndTracking.m and microrheology_1P.m have been run already
% and that the tracks for the individually tracked beads exist in a folder that is already
% given as a local variable
%
% Created by Daniel Seara at 2017/03/08 18:31
cd(basepath)

filelist = getAllFiles(pwd);

for ii = 1:numel(filelist)

end