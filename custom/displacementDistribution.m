% This script gives a histogram of the displacements at the tau that maximizes
% the non-Gaussian parameter found from nonGaussian.m at times before the cutoff
%
% alpha2(tau) = <dr^4>/((2) <dr^2>^2) - 1
%
% Assumes that featureFindingAndTracking.m and microrheology_1P.m have been run already
% and that the tracks for the individually tracked beads exist.
%
% INPUTS    basepath: path to msd folder after running all tracking and microrheology_1P scripts.
%                     Ends in filesep
%               tau : The time displacements used before the time cutoff to calculate the msd and non-
%                     Gaussian parameters, i.e. tau.pre
%            alpha2 : The non-Gaussian parameter before the time cutoff calculated with nonGaussian.m.
%                     i.e. alpha2.pre
%
% OUTPUTS         p : List of all displacements of all particles at the tau the maximizes alpha2
%
% Created by Daniel Seara at 2017/03/15 18:44

function p = displacementDistribution(basepath, tau, alpha2)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Preallocate memory and parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    preLength = length(tau);
    [~,maxFrame] = max(alpha2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ispc
        load([basepath 'Bead_Tracking\ddposum_files\individual_beads\correspondance_RG'])
    elseif isunix
        load([basepath 'Bead_Tracking/ddposum_files/individual_beads/correspondance_RG'])
    end

    p.x = [];
    p.y = [];

    for ii = 1:length(correspondance(:,1)) % Begin loop over beads

        if ispc
            load([basepath 'Bead_Tracking\ddposum_files\individual_beads\bead_' num2str(ii)]);
        elseif isunix
            load([basepath 'Bead_Tracking/ddposum_files/individual_beads/bead_' num2str(ii)]);
        end

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
    end % end loop over beads

%histogram(p,'Normalization','pdf')

