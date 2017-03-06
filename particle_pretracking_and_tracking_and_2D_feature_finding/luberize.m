
function [newtracks] = luberize(tracks)

% reassigns the unique ID# to 0,1,2,3...
% /presort will sort on ID# first, then reassign
% start will begin with that ID#

% function returns a new track array

ndat=length(tracks(1,:));
% if (keyword_set(presort)) then begin
%     newtracks=tracks(*,sort(tracks(ndat,*)))
% endif else begin
%     newtracks=tracks
% endelse

newtracks=tracks;

u=unq((newtracks(:,ndat))',[]);
ntracks=length(u);
u=[0,u];
for i=1:ntracks,  newtracks(u(i)+1:u(i+1),ndat) = i; end

% if (keyword_set(start)) then newtracks(ndat,*)=newtracks(ndat,*)+start

