function[] = mpretrack(basepath,fovn,featuresize,barrI,barrRg,barrCc,IdivRg,numframes,masscut,Imin,field)

% This program should be used when you have determined the values of its
% parameters using mpretrack_init. The calling sequence is
% essentially the same. The features and found by calling feature2D (with
% parameters other than feature size hardcoded).
%
% Note: feature2D requires the Image processing toolbox due to a call to
% imdilate in the localmax subfunction. Use feature2D_nodilate for an
% alternative which works almost as well.
%
% INPUTS :
% basepath - The base path for the experiment. "fov#_times.mat" files
%           should be there, and individual images should be in
%           "fov#\fov#_####.tif"
% fovn - ID# for the series of images (typically, one field of view)
% featuresize - The size of the feature you want to find.
% barrI - The minimum intensity you want to accept.
% barrRg - The maximum Rg squared you want to accept.
% barrCc - The maximum eccentricity you want to accept.
% IdivRg - minimum ratio of Intensity/pixel to be accepted (integrated
%           intensity / Rg squared of feature)
% numframes - The number of images you have in your series
% Imin - (optional) the minimum intensity for a pixel to be considered as a potential
%           feature.
% masscut - (optional) the masscut parameter for feature2D to remove false positives
%           before rifining the position to speed up the code.
% field - (optional) set to 0 or 1 if image is actually odd or even field of an interlaced 
%           image. All the masks will then be constructed with a 2:1 aspect ratio. Otherwise 
%           set to 2 for progressive scan cameras. Defaults to 2.
%
% Commented out:
% Inv - A logical for inverting the image (1 inverts, 0 doesn't) to look
%           for dark features instead of bright ones.
%
% Also, the program looks for the files "fov#_times.mat" for the "time"
% variable and the images files "fov#\fov#_####.tif" from the basepath.
%
% OUTPUTS
%
% - Creates a subfolder called "Feature_Finding" where it outputs the
% accepted features' MT matrix (from feature2D) as "MT_##_Feat_Size_"
% featuresize ".mat"
% - copies last frame into same subfolder
%
% MT matrix format:
% 1 row per bead per frame, sorted by frame number then x position (roughly)
% columns:
% 1:2 - X and Y positions (in pixels)
% 3   - Integrated intensity
% 4   - Rg squared of feature
% 5   - eccentricity
% 6   - frame #
% 7   - time of frame
%
% REVISION HISTORY
% written by Paul Fournier and Vincent Pelletier (Maria Kilfoil's group),
% latest revision 10/18/07
% 10/26/07 Vincent -- commented out the Inv keyword, added a ratio of
% Iint to Rg parameter
% 12/21/07 Maria -- added optional field

if nargin < 11, field = 2; end
if nargin < 10, Imin = 0; end
if nargin < 9, masscut = 0; end

tic,

pathin  = basepath;
if ispc
    pathout = [pathin 'Feature_finding\'];
elseif isunix
    pathout = [pathin 'Feature_finding/'];
end
[status, message, messageid] = mkdir( pathout );

d=0;
load([pathin 'fov' num2str(fovn) '_times.mat']);

for x = 1:numframes
    if ispc
        strnam=[pathin 'fov' num2str(fovn) '\fov' num2str(fovn) '_' num2str(x,'%04i') '.tif'];
    elseif isunix
        strnam=[pathin 'fov' num2str(fovn) '/fov' num2str(fovn) '_' num2str(x,'%04i') '.tif'];
    end
    img=imread(strnam);
%     if Inv == 1
%         img=255-img;
%     end

    M = feature2D(img,1,featuresize,masscut,Imin,field);

    if mod(x,50) == 0
        disp(['Frame ' num2str(x)])
        % partway save, useful if the computer tends to crash for some
        % reason
%         save([pathout 'MT_' num2str(fovn) '_Feat_Size_' num2str(featuresize) '_partial.mat'] ,'MT')
    end

    [a,b]=size(M);
    
    if( b == 5 )    

            %Rejection process
        X=find(M(:,5)>barrCc);
        M(X,1:5)=0;
        X=find(M(:,4)>barrRg);
        M(X,1:5)=0;
        X=find(M(:,3)<barrI);
        M(X,1:5)=0;
        X=find(M(:,3)./M(:,4)<IdivRg);
        M(X,1:5)=0;

        M=M(M(:,1)~=0,:);

        a = length(M(:,1));

        MT(d+1:a+d, 1:5)=M(1:a,1:5);
        MT(d+1:a+d, 6)=x;
        MT(d+1:a+d, 7)=time(x);
        d = length(MT(:,1));
        disp([num2str(a) ' features kept.'])
    end
    
    clear img;
    clear R;
    clear M;
    clear pic;
    clear X;
    clear t;
    clear i;
    clear j;
end

format long e;
% if Inv == 0
    save([pathout 'MT_' num2str(fovn) '_Feat_Size_' num2str(featuresize)],'MT')
% elseif Inv == 1
%     save([pathout 'MT_' num2str(fovn) '_Feat_Size_' num2str(featuresize) '_inv'],'MT')
% end
copyfile( strnam, [pathout 'fov' num2str(fovn) '_last.tif'] );    


clear all;
format short;

disp(['The program ran for ' num2str(toc/60) ' minutes'])
