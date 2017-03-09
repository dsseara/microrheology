function [M2,MT] = feat_find_init(basepath, featsize, barint, barrg, barcc, IdivRg, run, frame, masscut, Imin)
 
% This is an adapted version of the feature finding programs we have. It is
% only used to find the parameters (Intensity, Rg, Eccentricity, Feature
% Size). These parameters should then be plugged in
% feat_find_loop_p_fun_saving.
%
% Note: feature2D requires the Image processing toolbox due to a call to
% imdilate in the localmax subfunction. Use feature2D_nodilate for an
% alternative which works almost as well.
%
% INPUTS :
% basepath - The base path for the experiment. "run#_times.mat" files
% should be there, and individual images should be in
% "run#\run#_####.tif"
%
% featsize - The size of the feature you want to find.
% barint - The minimum intensity you want to accept.
% barrg - The maximum Rg squared you want to accept.
% barcc - The maximum eccentricity you want to accept.
% IdivRg - minimum ratio of Intensity/pixel to be accepted (integrated
% intensity / Rg squared of feature)
% run - ID# for the run
% frame - ID# for the image in the run
% Imin - (optional) the minimum intensity for a pixel to be considered as a potential
% feature (passed to feature2D)
% masscut - (optional) the masscut parameter for feature2D to remove false positives
% before rifining the position to speed up the code.
%
% commented out:
% inv - A logical for inverting the image (1 inverts, 0 doesn't) to look
% for dark features instead of bright ones.
%
% Also, the program looks for the files "run#_times.mat" for the "time"
% variable and the images files "run#\run#_####.tif" from the basepath.
%
% OUTPUTS 
% 
% M2 - All the features found from calling Nfeature3
% MT - All the features from Nfeature3 which were accepted given the
% criterions from the inputs
%
% The following appear on the command line:
%
% - Number of features found (without eliminating anything).
% - A matrix with the positions in columns 1-2 and the parameters int,rg
% and cc in columns 3-5. (accepted features first, then rejected features)
% - Min I, Max Rg, Max Cc these values help you sum up quickly what is in
% the matrix that appears. 
% - A plot showing the image and dots for rejected features and circles for
% accepted features
%
% REVISION HISTORY
% written by Paul Fournier and Vincent Pelletier (Maria Kilfoil's group),
% latest revision 10/18/07
% 10/26/07 Vincent -- commented out the Inv keyword, added a ratio of
% Iint to Rg parameter

if nargin < 10
    Imin = 0;
end
if nargin < 9
    masscut = 0;
end


pathin  = basepath;

d=0;
MT=[];
Mrej=[];

load([pathin 'run' num2str(run) '_times.mat']);
img=imread([pathin 'run' num2str(run) '\run' num2str(run) '_' num2str(frame,'%04i') '.tif']);

%     if inv == 1
%         img2=img;
%         img = 255-img;
%     end

    M = feature2D(img,1,featsize,masscut,Imin);

    a = length(M(:,1));
    
    M2=M;
     for i=1:a
         if ((M(i,5)>barcc))      
            Mrej=[Mrej; M(i,:)];
             M(i,1:5)=0;
%          end

         elseif ((M(i,4)>barrg))
            Mrej=[Mrej; M(i,:)];
             M(i,1:5)=0;
%          end

         elseif ((M(i,3)<barint))
            Mrej=[Mrej; M(i,:)];
             M(i,1:5)=0;
         elseif((M(i,3)/M(i,4)<IdivRg))
            Mrej=[Mrej; M(i,:)];
             M(i,1:5)=0;
         end    
     end
    
%    Deleting the zero rows

    M=M(M(:,1)~=0,:);

    a = length(M(:,1));

    MT(d+1:a+d, 1:5)=M(1:a,1:5);
%     MT(d+1:a+d, 6)=frame;
%     MT(d+1:a+d, 7)=time(frame);

    figure;
    imagesc(img,[0,255]),colormap(gray);
    hold on
    
    % Making a circle the size of the feature around each feature.
    
    theta = 0:0.001:2*pi;
    for c = 1:length(M(:,1))
        cx = M(c,1) + featsize*cos(theta)*2;
        cy = M(c,2) + featsize*sin(theta)*2;
        plot(cx,cy,'g-','linewidth',1.5)
    end
    if( ~isempty(Mrej)>0 )
        plot( Mrej(:,1), Mrej(:,2), 'r.' );
    end

%     if inv == 1
%         figure;
%         imagesc(img2,[0,255]),colormap(gray);
%         hold on
% 
%         % Making a circle the size of the feature around each feature.
%         for c = 1:length(M(:,1))
%             cx = M(c,1) + featsize*cos(theta);
%             cy = M(c,2) + featsize*sin(theta);
%             plot(cx,cy,'g-','linewidth',1.5)
%         end
%         if( ~isempty(Mrej)>0 )
%             plot( Mrej(:,1), Mrej(:,2), 'r.' );
%         end
%     end
    axis equal;
    
    format short g
    disp(M)
    disp(['Kept : ' num2str(size(M,1))])
    disp(Mrej)
    disp(['Minimum Intensity : ' num2str(min(M(:,3)))])
    disp(['Maximum Rg : ' num2str(max(M(:,4)))])
    disp(['Maximum Eccentricity : ' num2str(max(M(:,5)))])
    clear img;
    clear R;
    clear pic;
    clear X;
    clear t;
    clear i;
    clear j;



