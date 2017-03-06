function[] = conversions_no_dd(basepath, take)

% This program takes in the original res files (from trackmem) and converts 
% the positions in micrometers (no dedrifting, conversion factor 
% hardcoded in "pixtomicro"). (Follows "make_res_files")
%
% INPUTS
%
% basepath - The base path for the experiment. The function looks for the
% res matrices as "Bead_Tracking\res_files\res_run##.mat" files
% take - ID# for the run to be analyzed
%
% OUTPUTS
%
% creates subfolder "Bead_tracking\ddposum_files\" in which it saves
% "ddposum_run##.mat" files, with a truncated version of the "res" matrix
% The missing frames (possible when the "memory" parameter of trackmem is
% set to non-zero) are inserted as extrapolated positions.
%
% ddposum matrix format:
% 1 row per bead per frame, sorted by bead ID then frame number
% columns:
% 1:2 - X and Y positions (in micrometers)
% 3   - frame #
% 4   - Bead ID

if ispc
    pathout = ([basepath 'Bead_tracking\ddposum_files\']);
elseif isunix
    pathout = ([basepath 'Bead_tracking/ddposum_files/']);
end

[status,message,messageid] = mkdir( pathout );

if ispc
    load( [basepath 'Bead_Tracking\res_files\res_fov' num2str(take) '.mat'] );   % Load the tracked data as "res"
elseif isunix
    load( [basepath 'Bead_Tracking/res_files/res_fov' num2str(take) '.mat'] );   % Load the tracked data as "res"
end

res_acc = from_8_columnns_to_4(res);  % This part makes the matrix of 4 
                                       % columns that is more convenient.
                                       
ddposum = pixtomicro(res_acc);                     %converts it from pixels to microns
ddposum = putting_in_missing_frames(ddposum);

save( [pathout 'ddposum_run' num2str(take) '.mat'], 'ddposum' );
