function [MSD, tau, errMSD] = Mean_SD_many_single_beads_E(path, timeint, number_of_frames, rg_cutoff)
% This program calculates the mean squared displacement for all beads at 
% all possible lag times from 1 to number_of_frames. Also computes the
% standard error from the standard deviation at each tau.
% Follows "rg_matrix_many_single_beads".
%
% INPUTS
%
% path - the base path for the experiment. Reads in the individual
% beads files and the correspondance matrix with Rg from the 
% "Bead_Tracking\ddposum_files\individual_beads\" subfolder.
% timeint - the (average) delta-t between frames, in seconds
% number_of_frames - maximum lag time to consider (typically number of 
% frames recorded per FOV).
% rg_cutoff - maximum Rg value to be included. Beads with larger Rg are
% probably not confined by what we want to probe.
%
% OUTPUTS
%
% Creates a "1pt_msd" subfolder where it outputs
% "MSD_of_#_beads_rgcutoff_#nm.mat" and "MDSx..." and "MSDy..." files.
% (Usually not used)
% MSD - a vector containing the MSD in micrometer squared for each lag time tau
% tau - a vector of each lag time tau, in seconds
% errMSD - a vector containing the standard error in micrometer squared of 
% the MSD at each lag time tau
%

load([path 'Bead_Tracking\ddposum_files\individual_beads\correspondance_rg'])

msdpath = '1pt_msd';
[status, message, messageid] = mkdir([path msdpath]);

tau=timeint:timeint:(number_of_frames-1)*timeint;
MSD=zeros(number_of_frames-1,1);
MSDx=zeros(number_of_frames-1,1);
MSDy=zeros(number_of_frames-1,1);
sqMSD=zeros(number_of_frames-1,1);
sqMSDx=zeros(number_of_frames-1,1);
sqMSDy=zeros(number_of_frames-1,1);

pointtracer=zeros(number_of_frames-1,1);
beadcount = 0;

for i = 1:length(correspondance(:,1))
    clear M
    clear Mx
    clear My
    clear SD
    if correspondance(i,4) < rg_cutoff
        beadcount = beadcount + 1;
        load([path 'Bead_Tracking\ddposum_files\individual_beads\bead_' num2str(i)])
%         lastframe = max(bsec(:,3)); 
        lastframe=length(bsec(:,3));
        bsectauX=zeros(number_of_frames-1,1);
        bsectauY=zeros(number_of_frames-1,1);
        bsecx=(bsec(:,1)-bsec(1,1));
        bsecy=(bsec(:,2)-bsec(1,2));
   
        for delta=1:(lastframe-1)
            for k=1:(lastframe-delta)
                bsectauX(delta)=bsectauX(delta)+(bsecx(k)-bsecx(k+delta))^2;
                bsectauY(delta)=bsectauY(delta)+(bsecy(k)-bsecy(k+delta))^2;
                pointtracer(delta) = pointtracer(delta)+1; 
            end
        end

        MSDx=MSDx+bsectauX;
        MSDy=MSDy+bsectauY;
        MSD=MSD+bsectauX+bsectauY;
    end
    if mod(i,50) == 0
        disp(['Finished computing MSD for bead number ' num2str(i)])
    end

end


    
%howmanyPoints=numbeads*(number_of_frames-(1:number_of_frames-1))';

MSD=MSD./pointtracer;
MSDx=MSDx./pointtracer;
MSDy=MSDy./pointtracer;

for i = 1:length(correspondance(:,1))
    if correspondance(i,4) < rg_cutoff
        load([path 'Bead_Tracking\ddposum_files\individual_beads\bead_' num2str(i)])
%         lastframe = max(bsec(:,3)); 
        lastframe=length(bsec(:,3));
        sqbsectauX=zeros(number_of_frames-1,1);
        sqbsectauY=zeros(number_of_frames-1,1);
        bsecx=(bsec(:,1)-bsec(1,1));
        bsecy=(bsec(:,2)-bsec(1,2));    

        for delta=1:(lastframe-1)
            for k=1:(lastframe-delta)
                sqbsectauX(delta)=sqbsectauX(delta)+((bsecx(k)-bsecx(k+delta))^2-MSDx(delta))^2;
                sqbsectauY(delta)=sqbsectauY(delta)+((bsecy(k)-bsecy(k+delta))^2-MSDy(delta))^2;
                
            end
        end

        sqMSDx=sqMSDx+sqbsectauX;
        sqMSDy=sqMSDy+sqbsectauY;
        sqMSD=sqMSD+sqbsectauX+sqbsectauY;
    end
    if mod(i,50) == 0
        disp(['Finished computing sqMSD for bead number ' num2str(i)])
    end
end

sqMSD=sqrt(sqMSD./(pointtracer-1));     % standard deviation
sqMSDx=sqrt(sqMSDx./(pointtracer-1));
sqMSDy=sqrt(sqMSDy./(pointtracer-1));

errMSD=sqMSD./sqrt(pointtracer);    % standard error

save([path msdpath '\MSD_of_' num2str(beadcount) '_beads_rgcutoff_' num2str(rg_cutoff*1000) 'nm'],'MSD')
save([path msdpath '\MSDx_of_' num2str(beadcount) '_beads_rgcutoff_' num2str(rg_cutoff*1000) 'nm'],'MSDx')
save([path msdpath '\MSDy_of_' num2str(beadcount) '_beads_rgcutoff_' num2str(rg_cutoff*1000) 'nm'],'MSDy')
