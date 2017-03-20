% This script calculates a Non-gaussian parameter for all taus available
%
% alpha2(tau) = <dr^4>/((2) <dr^2>^2) - 1
%
% Assumes that featureFindingAndTracking.m and microrheology_1P.m have been run already
% and that the tracks for the individually tracked beads exist
%
% INPUTS     basepath       : path to data, ends in filesep
%            msd            : mean square displacement to use
%            tau            : time points to use
%            prepost_or_all : string, either 'prepost' or 'all' if going to use pre&post
%                              a critical time or just the entire series, respectively
%
% Created by Daniel Seara at 2017/03/08 18:31
function alpha2 = nonGaussian(basepath, msd, prepost_or_all)

    if ispc
        load([basepath 'Bead_Tracking\ddposum_files\individual_beads\correspondance_RG'])
    elseif isunix
        load([basepath 'Bead_Tracking/ddposum_files/individual_beads/correspondance_RG'])
    end
    
    switch prepost_or_all
        %%% Start pre-post analysis
        case 'prepost'
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Preallocate memory and parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            preLength = length(msd.pre);
            pre_pointtracer=zeros(preLength,1);
            mfd.pre=zeros(preLength,1); % Mean "fourth" displacement, <(delta x)^4>

            postLength = length(msd.post);
            post_pointtracer=zeros(postLength,1);
            mfd.post=zeros(postLength,1); % Mean "fourth" displacement, <(delta x)^4>
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            for ii = 1:length(correspondance(:,1)) % Begin loop over beads

                if ispc
                    load([basepath 'Bead_Tracking\ddposum_files\individual_beads\bead_' num2str(ii)]);
                elseif isunix
                    load([basepath 'Bead_Tracking/ddposum_files/individual_beads/bead_' num2str(ii)]);
                end

               %%% Pre tc %%%
                pre  = bsec(bsec(:,3)<(n+1),:);

                if isempty(pre)
                    % msdx.pre = 0;
                    % msdy.pre = 0;
                    % msd.pre  = 0;
                    disp('empty pre')
                    continue
                end

                pre_lastframe=length(pre(:,3));
                pre_bsectauX=zeros(preLength,1);
                pre_bsectauY=zeros(preLength,1);
                pre_bsecx=(pre(:,1)-pre(1,1));
                pre_bsecy=(pre(:,2)-pre(1,2));

                for delta=1:(pre_lastframe-1)
                    for k=1:(pre_lastframe-delta)
                        pre_bsectauX(delta) = pre_bsectauX(delta)+(pre_bsecx(k)-pre_bsecx(k+delta))^2;
                        pre_bsectauY(delta) = pre_bsectauY(delta)+(pre_bsecy(k)-pre_bsecy(k+delta))^2;
                        pre_pointtracer(delta) = pre_pointtracer(delta)+1; 
                    end
                end


                mfd.pre  = mfd.pre  + (pre_bsectauX + pre_bsectauY).^2;

                %%% Post tc %%%
                post  = bsec(bsec(:,3)>preLength,:);
                
                if isempty(post)
                    % msdx.post = 0;
                    % msdy.post = 0;
                    % msd.post  = 0;
                    disp('empty post')
                    continue
                end

                post_lastframe=length(post(:,3));
                post_bsectauX=zeros(postLength,1);
                post_bsectauY=zeros(postLength,1);
                post_bsecx=(post(:,1)-post(1,1));
                post_bsecy=(post(:,2)-post(1,2));

                for delta=1:(post_lastframe-1)
                    for k=1:(post_lastframe-delta)
                        post_bsectauX(delta)=post_bsectauX(delta)+(post_bsecx(k)-post_bsecx(k+delta))^2;
                        post_bsectauY(delta)=post_bsectauY(delta)+(post_bsecy(k)-post_bsecy(k+delta))^2;
                        post_pointtracer(delta) = post_pointtracer(delta)+1; 
                    end
                end

                mfd.post  = mfd.post  + (post_bsectauX + post_bsectauY).^2;
            end % end loop over beads
            
            mfd.pre  = mfd.pre ./ pre_pointtracer;
            mfd.post = mfd.post./ post_pointtracer;

            alpha2.pre  = (mfd.pre) ./((2) .* (msd.pre).^2)  - 1;
            alpha2.post = (mfd.post)./((2) .* (msd.post).^2) - 1;
        %%% end pre-post analysis

        %%% start all analysis
        case 'all'
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Preallocate memory and parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            n = length(msd);
            pointtracer=zeros(n,1);
            mfd=zeros(n,1); % Mean "fourth" displacement, <(delta x)^4>
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            for ii = 1:length(correspondance(:,1)) % Begin loop over beads

                if ispc
                    load([basepath 'Bead_Tracking\ddposum_files\individual_beads\bead_' num2str(ii)]);
                elseif isunix
                    load([basepath 'Bead_Tracking/ddposum_files/individual_beads/bead_' num2str(ii)]);
                end

                lastframe=length(bsec(:,3));
                bsectauX=zeros(n,1);
                bsectauY=zeros(n,1);
                bsecx=(bsec(:,1)-bsec(1,1));
                bsecy=(bsec(:,2)-bsec(1,2));

                for delta=1:(lastframe-1)
                    for k=1:(lastframe-delta)
                        bsectauX(delta) = bsectauX(delta)+(bsecx(k)-bsecx(k+delta))^2;
                        bsectauY(delta) = bsectauY(delta)+(bsecy(k)-bsecy(k+delta))^2;
                        pointtracer(delta) = pointtracer(delta)+1; 
                    end
                end


                mfd = mfd + (bsectauX + bsectauY).^2;
            end % end loop over beads
            
            mfd = mfd ./ pointtracer;

            alpha2 = mfd./((2) .* (msd).^2) - 1;
        %%% end all analysis
    end

