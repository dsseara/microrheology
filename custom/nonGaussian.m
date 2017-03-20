% This function calculates a Non-gaussian parameter for all taus available
% Assumes that featureFindingAndTracking.m and microrheology_1P.m have been run already
% and that the tracks for the individually tracked beads exist
%
% alpha2 = nonGaussian(basepath, msd, rg_cutoff, prepost_or_all) 
%
% INPUTS     basepath       : path to data, ends in filesep
%            msd            : Previously calculated mean square displacement
%            rg_cutoff      : [min max] of radius of gyration to use (same as for msd)
%            prepost_or_all : string, either 'prepost' or 'all' if going to use pre&post
%                              a critical time or just the entire series, respectively
%
% OUTPUTS    alpha2 : 2D non-Gaussian parameter given by  <dr^4>/((2) <dr^2>^2) - 1
%
% Created by Daniel Seara at 2017/03/20 01:54
function alpha2 = nonGaussian(basepath, msd, rg_cutoff, prepost_or_all)

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

                if correspondance(ii,4) < rg_cutoff(2) && correspondance(ii,4)>rg_cutoff(1)
                   %%% Pre tc %%%
                    pre  = bsec(bsec(:,3)<(preLength+1),:);

                    if isempty(pre)
                        % msdx.pre = 0;
                        % msdy.pre = 0;
                        % msd.pre  = 0;
                        disp('empty pre')
                        continue
                    end

                    pre_lastframe=length(pre(:,3));
                    pre_bsectauR=zeros(preLength,1);
                    pre_bsecx=(pre(:,1)-pre(1,1));
                    pre_bsecy=(pre(:,2)-pre(1,2));

                    for delta=1:(pre_lastframe-1)
                        for k=1:(pre_lastframe-delta)
                            pre_bsectauR(delta) = pre_bsectauR(delta) + ((pre_bsecx(k)-pre_bsecx(k+delta))^2 + (pre_bsecy(k)-pre_bsecy(k+delta))^2)^2;
                            pre_pointtracer(delta) = pre_pointtracer(delta)+1; 
                        end
                    end


                    mfd.pre  = mfd.pre  + pre_bsectauR;

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
                    post_bsectauR=zeros(postLength,1);
                    post_bsecx=(post(:,1)-post(1,1));
                    post_bsecy=(post(:,2)-post(1,2));

                    for delta=1:(post_lastframe-1)
                        for k=1:(post_lastframe-delta)
                            post_bsectauR(delta) = post_bsectauR(delta) + ((post_bsecx(k)-post_bsecx(k+delta))^2 + (post_bsecy(k)-post_bsecy(k+delta))^2)^2;
                            post_pointtracer(delta) = post_pointtracer(delta)+1; 
                        end
                    end

                    mfd.post  = mfd.post  + post_bsectauR;
                end
            end % end loop over beads
            
            mfd.pre  = mfd.pre ./ pre_pointtracer;
            mfd.post = mfd.post./ post_pointtracer;

            alpha2.pre  = (mfd.pre) ./(2.*(msd.pre).^2)  - 1;
            alpha2.post = (mfd.post)./(2.*(msd.post).^2) - 1;
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
                
                if correspondance(ii,4) < rg_cutoff(2) && correspondance(ii,4)>rg_cutoff(1)
                    
                    lastframe=length(bsec(:,3));
                    bsectauR=zeros(n,1);
                    bsecx=(bsec(:,1)-bsec(1,1));
                    bsecy=(bsec(:,2)-bsec(1,2));

                    for delta=1:(lastframe-1)
                        for k=1:(lastframe-delta)
                            bsectauR(delta) = bsectauR(delta) + ((bsecx(k)-bsecx(k+delta))^2 + (bsecy(k)-bsecy(k+delta))^2)^2;
                            pointtracer(delta) = pointtracer(delta)+1; 
                        end
                    end

                    mfd = mfd + bsectauR;
                end
            end % end loop over beads
            
            mfd = mfd ./ pointtracer;
            alpha2 = mfd./(2.*(msd).^2) - 1;
        %%% end all analysis
    end

