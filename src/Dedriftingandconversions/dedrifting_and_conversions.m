function[] = dedrifting_and_conversions(basepath, take, smoo)

% This program takes in the original res files (from trackmem), removes its 
% drift and converts the positions in micrometers (conversion factor 
% hardcoded in "pixtomicro"). (Follows "make_res_files")
%
% INPUTS
%
% basepath - The base path for the experiment. The function looks for the
% res matrices as "Bead_Tracking\res_files\res_run##.mat" files
% take - ID# for the run to be analyzed
% smoo - the number of frames over which the drift is averaged before being
% removed from the bead positions.
%
% OUTPUTS
%
% creates subfolder "Bead_tracking\ddposum_files\" in which it saves
% "ddposum_run##.mat" files, with a truncated version of the "res" matrix
% The missing frames (possible when the "memory" parameter of trackmem is
% set to non-zero) are inserted as extrapolated positions.
% Also saves a plot of the drift and dedrift line as "dedrift_run##.png"
% (time axis hardcoded in "dirft_loop").
%
% ddposum matrix format:
% 1 row per bead per frame, sorted by bead ID then frame number
% columns:
% 1:2 - X and Y positions (in micrometers)
% 3   - frame #
% 4   - Bead ID

pathout = ([basepath 'Bead_tracking\ddposum_files\']);   
[status,message,messageid] = mkdir( pathout );

load( [basepath 'Bead_Tracking\res_files\res_run' num2str(take) '.mat'] );   % Load the tracked data as "res"
ddposum = pixtomicro(res);      %converts it from pixels to microns, using hardcoded values

ddposum = drift_loop(ddposum,smoo);                      %finds and removes the drift
print('-dpng','-r150', [pathout 'dedrift_run' num2str(take) '.png'])
close all


save( [pathout 'ddposum_run' num2str(take) '.mat'], 'ddposum' );
