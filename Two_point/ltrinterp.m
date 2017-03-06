function []=ltrinterp(tarray)

%%%%get back to this function later

% ;   started 6-17-98 by ERW, modified for use here 6/99 by JCC
% ;   interpolates gaps in tracked data (time,id only), zeros elsewhere
% ;   tarray is an array of tracked data

ndat=length(tarray(1,:));
dp=tarray-circshift(tarray,[1,0]);
w=find((dp(:,ndat) == 0) && (abs(dp(:,ndat-1)-1) > 1e6); ngood=length(w);
if (ngood > 1)
    totdt=sum(dp(w,ndat-1))-ngood;
    dp=0;
    storeres=zeros(totdt,ndat);
    count = 0;
    for i=1:ngood
       dt = tarray(w(i),ndat-1) - tarray(w(i),ndat-1);
       timer = tarray(w(i),ndat-1) + findgen(dt,2) +1;
       storeres(ndat-1,count:count+dt-2) = tarray(ndat-1,w(i)-1)
       storeres(ndat-2,count:count+dt-2) = timer
       count = count + dt - 1
    end
    tarray = [[tarray],[storeres]]
    ; watch out on other platforms' sorts!
    tarray = tarray(*,sort(tarray(ndat-2,*)))
    tarray = tarray(*,sort(tarray(ndat-1,*)))
endif
end
;