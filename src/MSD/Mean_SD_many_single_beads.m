function [msd, msdx, msdy, tau, beadcount] = Mean_SD_many_single_beads(path, timeint, number_of_frames,rg_cutoff, tc)

% This program calculates the mean squared displacement for all beads at 
% all possible lag times from 1 to number_of_frames. The result is the sum
% of the x and y msd.
% Follows "rg_matrix_many_single_beads".
%
% INPUTS
%
%             path - the base path for the experiment. Reads in the individual
%                    beads files and the correspondance matrix with Rg from the 
%                    "Bead_Tracking\ddposum_files\individual_beads\" subfolder.
%          timeint - the (average) delta-t between frames, in seconds
% number_of_frames - maximum lag time to consider (typically number of 
%                    frames recorded per FOV).
%               tc - (optional) critical time (as a frame number) to calculate msd
%                    both before and after is set 
%
% OUTPUTS
%
% Creates a "1pt_msd" subfolder where it outputs
% "msd_of_#_beads_rgcutoff_#nm.mat" and "MDSx..." and "msdy..." files.
% (Usually not used)
% msd - a vector containing the msd in micrometer squared for each lag time tau
% tau - a vector of each lag time tau, in seconds
%

%load([path 'Bead_Tracking\ddposum_files\individual_beads\correspondance'])
if ispc
    load([path 'Bead_Tracking\ddposum_files\individual_beads\correspondance_RG'])
elseif isunix
    load([path 'Bead_Tracking/ddposum_files/individual_beads/correspondance_RG'])
end

if isempty(tc)
    tau=[timeint:timeint:(number_of_frames-1)*timeint]';
    msd=zeros(number_of_frames-1,1);
    msdx=zeros(number_of_frames-1,1);
    msdy=zeros(number_of_frames-1,1);
    pointtracer=zeros(number_of_frames-1,1);
else
    %%% Pre tc %%%
    pre_number_of_frames = tc;
    tau.pre=[timeint:timeint:(pre_number_of_frames-1)*timeint]';
    pre_pointtracer=zeros(pre_number_of_frames-1,1);
    %%% Post tc %%%
    post_number_of_frames = number_of_frames - tc;
    tau.post=[timeint:timeint:(post_number_of_frames-1)*timeint]';
    post_pointtracer=zeros(post_number_of_frames-1,1);

    % Make two structs, one for pre and one for post tc
    msd.pre=zeros(pre_number_of_frames-1,1);
    msdx.pre=zeros(pre_number_of_frames-1,1);
    msdy.pre=zeros(pre_number_of_frames-1,1);

    msd.post=zeros(post_number_of_frames-1,1);
    msdx.post=zeros(post_number_of_frames-1,1);
    msdy.post=zeros(post_number_of_frames-1,1);
end

if isempty(rg_cutoff)
    rg_cutoff = max(correspondance(:,4)) + 1;
end

beadcount = 0;

for i = 1:length(correspondance(:,1))
    clear M
    clear Mx
    clear My
    clear SD
    beadcount = beadcount + 1;
    
    if ispc
        load([path 'Bead_Tracking\ddposum_files\individual_beads\bead_' num2str(i)]);
    elseif isunix
        load([path 'Bead_Tracking/ddposum_files/individual_beads/bead_' num2str(i)]);
    end
    
    if correspondance(i,4) < rg_cutoff
        if isempty(tc)
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

            msdx=msdx+bsectauX;
            msdy=msdy+bsectauY;
            msd=msd+bsectauX+bsectauY;
        else
            %%% Pre tc %%%
            pre  = bsec(bsec(:,3)<tc+1,:);
            
            if isempty(pre)
                % msdx.pre = 0;
                % msdy.pre = 0;
                % msd.pre  = 0;
                disp('empty pre')
                continue
            end

            pre_lastframe=length(pre(:,3));
            pre_bsectauX=zeros(pre_number_of_frames-1,1);
            pre_bsectauY=zeros(pre_number_of_frames-1,1);
            pre_bsecx=(pre(:,1)-pre(1,1));
            pre_bsecy=(pre(:,2)-pre(1,2));

            for delta=1:(pre_lastframe-1)
                for k=1:(pre_lastframe-delta)
                    pre_bsectauX(delta) = pre_bsectauX(delta)+(pre_bsecx(k)-pre_bsecx(k+delta))^2;
                    pre_bsectauY(delta) = pre_bsectauY(delta)+(pre_bsecy(k)-pre_bsecy(k+delta))^2;
                    pre_pointtracer(delta) = pre_pointtracer(delta)+1; 
                end
            end

            msdx.pre = msdx.pre + pre_bsectauX;
            msdy.pre = msdy.pre + pre_bsectauY;
            msd.pre  = msd.pre  + pre_bsectauX + pre_bsectauY;

            %%% Post tc %%%
            post  = bsec(bsec(:,3)>tc,:);
            
            if isempty(post)
                % msdx.post = 0;
                % msdy.post = 0;
                % msd.post  = 0;
                disp('empty post')
                continue
            end

            post_lastframe=length(post(:,3));
            post_bsectauX=zeros(post_number_of_frames-1,1);
            post_bsectauY=zeros(post_number_of_frames-1,1);
            post_bsecx=(post(:,1)-post(1,1));
            post_bsecy=(post(:,2)-post(1,2));

            for delta=1:(post_lastframe-1)
                for k=1:(post_lastframe-delta)
                    post_bsectauX(delta)=post_bsectauX(delta)+(post_bsecx(k)-post_bsecx(k+delta))^2;
                    post_bsectauY(delta)=post_bsectauY(delta)+(post_bsecy(k)-post_bsecy(k+delta))^2;
                    post_pointtracer(delta) = post_pointtracer(delta)+1; 
                end
            end

            msdx.post = msdx.post + post_bsectauX;
            msdy.post = msdy.post + post_bsectauY;
            msd.post  = msd.post  + post_bsectauX + post_bsectauY;
        end
    end
    if mod(i,50) == 0
        disp(['Finished computing msd for bead number ' num2str(i)])
    end
end


%howmanyPoints=numbeads*(number_of_frames-(1:number_of_frames-1))';
if isempty(tc)
    msd=msd./pointtracer;
    msdx=msdx./pointtracer;
    msdy=msdy./pointtracer;

    % save([saveTo 'wholeSeries_output'],'msd', 'msdx', 'msdy', 'tau')
else

    %%% Pre tc in first field %%%
    msd.pre =msd.pre ./pre_pointtracer;
    msdx.pre=msdx.pre./pre_pointtracer;
    msdy.pre=msdy.pre./pre_pointtracer;

    %%% Post tc in second field %%%
    msd.post =msd.post ./post_pointtracer;
    msdx.post=msdx.post./post_pointtracer;
    msdy.post=msdy.post./post_pointtracer;

    % save([saveTo 'tcFrame_' num2str(tc) '_output'], 'msd','msdx', 'msdy', 'tau')
end
