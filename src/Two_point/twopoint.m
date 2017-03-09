% ;function msd2pnt,tracks,params,[outfile=outfile,micperpix=micperpix,$
% ;    timestep=timestep,dim=dim,mydts=mydts,erode=erode,nodedrift=nodedrift]
% ;
% ;    PURPOSE:
% ;
% ;   MSD2PNT -- a 2 or 3 dimensional space-time software correlator
% ;
% ;   This is a little routine for measuring the 2-point mean
% ;   squared displacement from tracers in a viscoelastic medium.
% ;   This can then be plugged into 'micrheo2pnt' to get the medium's
% ;   G',G"(w).  It is strictly defined as <x1*x2>(dr,dt), where the
% ;   average is over time and ensembles.  The result triple counts
% ;   in time, is logarithmically spaced in 'dt' and 'dr', and can be
% ;   process 2 or 3 dimensional 'track'ed data.
% ;
%      INPUTS
%
%      tracks - the data matrix, containing positions, frame number and 
%      bead number (e.g. created by "separating_beads_by_FOV")
%      params - [rmin,rmax,nrbins,maxtime], determining the inner and
%      outer cutoffs for the radius, the number of log-spaced bins in r
%      and the maximum log-spaced time in seconds.  Lengths are in um.
%      timestep - number of seconds between frames
%      dim - spatial dimensionality of the data (2 or 3, though dim==3 has 
%      not been translated in that code yet) 
%      mydts - vector of dt's for which the correlations should be computed. 
%      Set to 0 to let the program create a logarithmically spaced vector.
%      dedrift - Set to 1 to perform dedrift, but really should be left to
%      0, we are doing the dedrift elsewhere already.
%
% ;    OUPUTS:
% ;   The result 'data' has the form:
% ;   In 2 dimensions (polar coordinates):
% ;     data(*,*,0) contains the log-spaced 'dr' values.
% ;     data(*,*,1) contains the log-spaced 'dt' values.
% ;     data(*,*,2) is the longitudinal, 'r-r' mean component.
% ;     data(*,*,3) is the transverse, 'theta-theta' mean component.
% ;     data(*,*,4:5) store the *variances* of the r,theta parts.
% ;     data(*,*,6) contains the total number of points used in each bin.
% ;   In 3 dimensions (spherical coordinates):
% ;     data(*,*,0) contains the log-spaced 'dr' values.
% ;     data(*,*,1) contains the log-spaced 'dt' values.
% ;     data(*,*,2) is the longitudinal, 'r-r' mean component.
% ;     data(*,*,3) is the transverse, 'theta-theta' mean component.
% ;     data(*,*,4) is the transverse, 'phi-phi' mean component.
% ;     data(*,*,5:7) store the *variances* of the r,theta,phi parts.
% ;     data(*,*,8) contains the total number of points used in each bin.
% ;
% ;    NOTES:
% ;     The means are what you want, the variances and numbers are
% ;     intended for use in error estimation of the mean values.
% ;     The means correspond to the diagonal components of a
% ;     correlation tensor that would be called 'S' in the
% ;     Fluctuation/Dissipation Theorem in Chaikin & Lubensky.
% ;     The units of 'S' are um^2 if the inputs are in um.
% ;
% ;    MODIFICATION HISTORY:
% ;         Written by John C. Crocker, U. Penn:    6/'99
% ;
% ;       This code 'msd2pnt.pro' is copyright 1999 by John C. Crocker
% ;       and may be freely distributed and used if properly attributed.

% ;   The actual routine itself.
% ;   'params' is [rmin,rmax,nrbins,dtmax].

function [data]=twopoint(tracks,params,timestep,dim,mydts,dedrift)
%more inputs can be added later: outfile, erode

%%%skipping this part, which allows tracks to be a file
% ; read the gdf file if 'tracks' is a string
% sz = size(tracks)
% nz = n_elements(sz)
% if sz(nz-2) eq 7 then begin
%     f = findfile(tracks)
%     if f(0) eq '' then message,'No Match: '+tracks
%     nf = n_elements(f)
%     filebased = 1
% endif else begin
%     filebased = 0
%     nf = 1
      erwmax=max(tracks(:,length(tracks(1,:))-1));
      erwmin=min(tracks(:,length(tracks(1,:))-1));
      params(4) = min( params(4) , (erwmax-erwmin)*timestep );
% endelse
% 

if mydts ==0
    % generate the time partition-- about 10 points per decade
    dt = round(1.15.^(0:99));
    [eldt,idt]=unique(dt);
    dt = dt(idt);
    round(params(4)/timestep);
    w = find( dt <= round(params(4)/timestep) ); ndt=length(w);
    if ndt > 0, dt = dt(w); else warning('Invalid maximum dt!'), end
else
    dt = mydts;
    ndt = length(dt);
end

%; generate the 'r' partition-- as the user wishes.
nbins = params(3);
lmax = log(params(2)); lmin = log(params(1));
lrbinsize = (lmax-lmin)/nbins;
rpart = exp( lmin+((0:nbins)*lrbinsize ));  % the r partitions
r = sqrt( rpart .* circshift(rpart,[0,-1]) );
r = r(1:nbins);            % the r bin 'midpoints'

%; set up some arrays-- the running sums of the data
data = zeros(ndt,nbins,(2*dim)+3);    %sL,sT,sL^2,sT^2,N(r,t) and t,r
%note: rows correspond to time, columns to distance. should be like in IDL

%; do the big nested loop
% for i=0,nf-1 do begin
%     if filebased then begin
%        message,'Reading file: '+f(i),/inf
%        trj = read_gdf(f(i))
%     endif else trj = tracks
%     ; convert to um.
%     for j=0,dim-1 do trj(j,*) = trj(j,*)*micperpix(j)
trj = tracks;
%%% trj=ltrinterp(trj);            % fill in gaps.   %%%include this later
for j=1:ndt
%        if filebased then $
%           message,'Processing file: '+f(i)+$
%          '    dt equals: '+strcompress(string(dt(j))),/inf $
%           else message,' dt equals: '$
%          +strcompress(string(dt(j))),/inf
    disp(['dt equals: ' num2str(dt(j)) ' j=' num2str(j)])
    res = get_corr(trj,[lmin,lmax,lrbinsize,nbins],dt(j),dim,dedrift);
%   ; accumulate the running sums
    res1(1,:,:)=res';
    data(j,:,1:(2*dim+1)) = data(j,:,1:(2*dim+1)) + res1;
    clear res1
end
% end

%from here on check again

%; make the running sums into means and variances.
for i = 1:(2*dim), data(:,:,i) = data(:,:,i)./data(:,:,2*dim+1); end
data(:,:,dim+1:(2*dim)) = data(:,:,dim+1:(2*dim)) - (data(:,:,1:dim).^2);

%; put on dt,dr vectors (slightly redundant data structure, but hey)
for i = 1:ndt, data(i,:,(2*dim)+2) = r; end
for i = 1:nbins, data(:,i,(2*dim)+3) = dt*timestep; end

%; shift the array so dr,dt is on the front like with msd.pro
data = circshift(data,[0,0,2]);

% %; optionally write out a file
% if keyword_set(outfile) then begin
%     message,'Writing file: '+outfile,/inf
%     write_gdf,data,outfile
% endif
