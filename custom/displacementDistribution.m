% This function gives histograms of the displacements at the tau that maximizes
% the non-Gaussian parameter found from nonGaussian.m at times before the cutoff
%
% alpha2(tau) = <dr^4>/((2) <dr^2>^2) - 1
%
% Assumes that featureFindingAndTracking.m and microrheology_1P.m have been run already
% and that the tracks for the individually tracked beads exist.
%
% p = displacementDistribution(basepath, tau, alpha2, rg_cutoff)
%
% INPUTS    basepath: path to msd folder after running all tracking and microrheology_1P scripts.
%                     Ends in filesep
%               tau : The time displacements used before the time cutoff to calculate the msd and non-
%                     Gaussian parameters, i.e. tau.pre
%            alpha2 : The non-Gaussian parameter before the time cutoff calculated with nonGaussian.m.
%                     i.e. alpha2.pre
%         rg_cutoff : [min max] of radius of gyration to use (same as for msd)
%
% OUTPUTS         p : Struct of displacements, cartesian (x, y, and combined), radial, and angular
%
% Created by Daniel Seara at 2017/03/15 18:44

function p = displacementDistribution(basepath, tau, alpha2, rg_cutoff)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Preallocate memory and parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    preLength = length(tau);
    [~,maxFrame] = max(alpha2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ispc
        load([basepath 'Bead_Tracking\ddposum_files\individual_beads\correspondance_rg.mat'])
    elseif isunix
        load([basepath 'Bead_Tracking/ddposum_files/individual_beads/correspondance_rg.mat'])
    end

    p.x = [];
    p.y = [];

    for ii = 1:length(correspondance(:,1)) % Begin loop over beads

        if ispc
            load([basepath 'Bead_Tracking\ddposum_files\individual_beads\bead_' num2str(ii)]);
        elseif isunix
            load([basepath 'Bead_Tracking/ddposum_files/individual_beads/bead_' num2str(ii)]);
        end

        if correspondance(ii,4) < rg_cutoff(2) && correspondance(ii,4)>rg_cutoff(1)
            %%% Pre tc %%%
            pre  = bsec(bsec(:,3)<(preLength+1),:);

            if isempty(pre)
                disp('empty pre')
                continue
            end

            pre_lastframe=length(pre(:,3));
            pre_bsectauX=zeros((pre_lastframe-maxFrame),1);
            pre_bsectauY=zeros((pre_lastframe-maxFrame),1);
            pre_bsecx=(pre(:,1)-pre(1,1));
            pre_bsecy=(pre(:,2)-pre(1,2));

            for k=1:(pre_lastframe-maxFrame)
                pre_bsectauX(k) = pre_bsectauX(k)+(pre_bsecx(k)-pre_bsecx(k+maxFrame));
                pre_bsectauY(k) = pre_bsectauY(k)+(pre_bsecy(k)-pre_bsecy(k+maxFrame));
            end


            %pre_bsectauR = (pre_bsectauX + pre_bsectauY).^(0.5);

            p.x = [p.x; pre_bsectauX];
            p.y = [p.y; pre_bsectauY];
        end
    end % end loop over beads

p.all = [p.x;p.y];
p.r = sqrt((p.x).^2 + (p.y).^2);
p.theta = atan2(p.y,p.x);

save('dispData.mat', 'p')




