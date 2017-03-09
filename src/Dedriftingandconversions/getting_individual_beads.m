function[] = getting_individual_beads( basepath, FOVs )

% This function is used to get position matrices for beads individually.
% This will simplify the use of later programs, such as seperating beads by
% rg and for finding MSD. Follows "conversions_no_dd_vp" or
% "dedrifting_and_conversions_large_smoo_dd".
%
% INPUTS
%
% basepath - the base path for the experiment. Reads in the ddposum matrix files for
% each field of view.
% FOVs - a vector of run ID # to be processed
%
% OUTPUTS
%
% creates a subfolder "Bead_tracking\ddposum_files\individual_beads\" in
% which the position matrix files are saved as "bead_###.mat".
% Saves a correspondance matrix which associates each bead # with it's FOV
% and bead ID in that FOV
%

if ispc
    pathin  = [ basepath '\Bead_tracking\ddposum_files\'];
    pathout = [pathin '\individual_beads\'];
elseif isunix
    pathin  = [ basepath '/Bead_tracking/ddposum_files/'];
    pathout = [pathin '/individual_beads/'];
end

[status,message,messageid] = mkdir( pathout );

h = 1;
for FOV = FOVs
    load([pathin 'ddposum_run' num2str(FOV) '.mat']);     %load ddposum
    poscor=ddposum;
    clear ddposum;

    beadindx = unique(fix(poscor(:,4)));
    beadmax=length(beadindx);
    for bead = 1:beadmax
        bsec = poscor(find(poscor(:,4)==beadindx(bead)),:);

        framemin = unique(min(bsec(:,3)));
        bsec(:,5) = bsec(:,3)- framemin + 1;        % why do we even need that?

        framemax = unique(max(bsec(:,3)));
        save( [pathout 'bead_' num2str(h) '.mat'],'bsec' );
        correspondanceh(h,1) = FOV;
        correspondanceh(h,2) = beadindx(bead);
        correspondanceh(h,3) = h;
        h = h+1;
    end
  disp(['Finished FOV ' num2str(FOV) '...'])
end

dlmwrite([pathout 'correspondance_3000.txt'],correspondanceh,'delimiter','\t','newline','pc')
correspondance=correspondanceh;
if ispc
    save([pathout '\correspondance.mat'],'correspondance');
elseif isunix
    save([pathout '/correspondance.mat'],'correspondance');
end


