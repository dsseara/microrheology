function result=getdx(tr,dt,inputv, d)
%program written by Yongxiang Gao on Nov 28,2005
%result contains the displacement vector for x y z and the |r|.
%It excludes those postion involves the interpolation.
if inputv==0; d=3;end
trin=gdxtrinterp(tr, 1,1); %what is flag? 
nel=length(trin(1,:));
trins=circshift(trin, [-dt,0]);
fl1 = trin(:,nel-2);
fl2 = trins(:,nel-2);
result=trins-trin;
w=find(result(:,nel)~=0 | result(:,nel-1)~=dt | fl1>0.5 | fl2>0.5);
nw=length(w);
result=trins(:,1:d+1)-trin(:,1:d+1);
trins=0; %free up memory
result(:,d+1)=sum(result(:,1:d).^2,2);
result(:,d+1)=sqrt(result(:,d+1));
if nw>0; result(w,:)=-1; end
w2=find(trin(:,nel-2) <0.5);
result=result(w2,:);
w3=find(result(:,d+1)<0);
nw3=length(w3);
if nw3>0; result(w3,1:d)=0; end
