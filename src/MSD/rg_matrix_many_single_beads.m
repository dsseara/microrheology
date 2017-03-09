function[] = rg_matrix_many_single_beads(basepath)

% This program gets the radius of gyration of all the different beads. 
% Follows "getting_individual_beads".
%
% INPUTS
%
% basepath - the base path for the experiment. Reads in the individual
% beads files and the correspondance matrix from the 
% "Bead_Tracking\ddposum_files\individual_beads\" subfolder.
%
% OUTPUTS
%
% in the "Bead_Tracking\ddposum_files\individual_beads\" subfolder, creates
% a new file "correspondance_rg.mat" which contains the correspondance
% matrix with the computed Rg.
%
% correspondance matrix format:
% 1 row per bead
% columns:
% 1 - FOV number
% 2 - Bead ID# in that FOV
% 3 - global bead #
% 4 - Rg in micrometers
%

%basepath = 'C:\Users\student\Kevin\Data_actin_patches\071213_actin_patches_KS\fov7\Analysis_0703_1_Cell1';
%load([basepath 'tracks_um_Vincent'])
if ispc
    load([basepath 'Bead_Tracking\ddposum_files\individual_beads\correspondance.mat']);
elseif isunix
    load([basepath 'Bead_Tracking/ddposum_files/individual_beads/correspondance.mat']);
end

for bead = 1:length(correspondance(:,1))
    if ispc
        load([basepath 'Bead_Tracking\ddposum_files\individual_beads\bead_' num2str(bead)])
    elseif isunix
        load([basepath 'Bead_Tracking/ddposum_files/individual_beads/bead_' num2str(bead)])
    end
    
    Rg = calc_rg_k(bsec,2);                                 %finds rg for the kth section 
    correspondance(bead,4) = Rg;
    clear Rg bsec
    
    if mod(bead,50) == 0
        disp(['Finished computing Rg for bead number ' num2str(bead)])
    end
end

if ispc
    save([basepath 'Bead_Tracking\ddposum_files\individual_beads\correspondance_rg'],'correspondance')
elseif isunix
    save([basepath 'Bead_Tracking/ddposum_files/individual_beads/correspondance_rg'],'correspondance')
end
