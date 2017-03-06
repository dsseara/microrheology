function result=gdxtrinterp(tarray, inputv,flag) %need to pay some attention on the flag

s=size(tarray);
result=tarray;
if inputv==1
    result=[result, zeros(s(1),1)];
end

ndat=length(tarray(1,:));
dp=tarray-circshift(tarray,[1 0]);
w=find(dp(:,ndat)==0 & dp(:,ndat-1)>1); %if w is empty, there is nothing we need to trinterp
ngood=length(w);
count=0;
storeres=0;
if ngood>=1
    totdt=sum(dp(w,ndat-1))-ngood;
    dp=0; %free up memory
    storeres=zeros(totdt, ndat);
    for i=1:ngood
        x0=tarray(w(i)-1,:);
        x1=tarray(w(i),:);
       dt=x1(ndat-1)-x0(ndat-1);
       t=1.0-(1:dt-1)/dt;
       xx = x0'*t+x1'*(1.0-t);   
       storeres(count+1:count+dt-1,:)=xx';
       count = count + dt-1;
    end
    n=length(storeres(:,1));
    if inputv==1
       storeres=[storeres,ones(n,1)];
    result=[[result];[storeres]];
    storeres=0;
    end
else
   ; %nothing to trinterp, apparently
end
% ; put in Victor Breedveld's patch here, to make sort work in DOS
[Y,id]=sort(result(:,ndat));
result=result(id,:);
[B, idx]=unique(result(:,ndat));
u=[-1; idx];
nu=length(u);
for i=2:nu
    tpart=result(u(i-1)+2:u(i),:);  
    [Y,id]=sort(tpart(:,ndat-1));
    result(u(i-1)+2:u(i),:)=tpart(id,:);
end
% ; end of Victor's patch
if (inputv==1)
    result(:,ndat-1:end)=result(:,[ndat+1, ndat-1, ndat]);
end
% ; this ends the function

