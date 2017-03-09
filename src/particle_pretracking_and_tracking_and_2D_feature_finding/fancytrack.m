function [ ] = fancytrack( basepath, FOVnum, featsize, maxdisp, goodenough, memory )

% Runs trackmem on the output from mpretrack. 
%
% INPUTS
%
% basepath - the basepath of the experiments. Reads the MT matrix from
%       "Feature_finding\MT_##_Feat_Size_" featsize ".mat", as output by
%       mpretrack
% FOVnum - specifies which series of images to process
% featsize - specifies the feature size for accessing the right MT file
% maxdisp - (optional) specifies the maximum displacement (in pixels) a feature may 
%       make between successive frames
% goodenough - (optional) the minimum length requirement for a trajectory to be retained 
% memory - (optional) specifies how many consecutive frames a feature is allowed to skip. 
%
% OUTPUTS
%
% creates subfolders "Bead_tracking\res_files\" in which it saves
% "res_fov##.mat" files, with a "res" matrix which is the output from
% trackmem
%
% res matrix format:
% 1 row per bead per frame, sorted by bead ID then by frame number.
% columns are:
% 1:2 - X and Y positions (in pixels)
% 3   - Integrated intensity
% 4   - Rg squared of feature
% 5   - eccentricity
% 6   - frame #
% 7   - time of frame
% 8   - Bead ID
%
% REVISION HISTORY
% written by Paul Fournier and Vincent Pelletier (Maria Kilfoil's group),
% last revision 10/18/07

if nargin < 6, memory=1; end
if nargin < 5, goodenough=100; end
if nargin < 4, maxdisp=2; end

[status,message,messageid] = mkdir( [basepath 'Bead_tracking'], 'res_files' );

j=FOVnum;

%     if( bUseInverted==1 )
%         load([basepath 'Feature_finding\MT_' num2str(j) '_Feat_Size_' num2str(featsize) '_inv.mat'])
%         MTinv=MT;
%     end
    if ispc
        load([basepath 'Feature_finding\MT_' num2str(j) '_Feat_Size_' num2str(featsize) '.mat'])
    elseif isunix
        load([basepath 'Feature_finding/MT_' num2str(j) '_Feat_Size_' num2str(featsize) '.mat'])
    end
    
%     if( bUseInverted==1 )
%         MT=[MT;MTinv];
%         MT=sortrows( MT, 6 );
%     end
    
    res=trackmem( MT, maxdisp, 2, goodenough, memory );
    if ispc
        save( [basepath 'Bead_tracking\res_files\res_fov' num2str(j) '.mat'], 'res' );
    elseif isunix
        save( [basepath 'Bead_tracking/res_files/res_fov' num2str(j) '.mat'], 'res' );
    end


