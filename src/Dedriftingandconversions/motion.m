function result=motion(tr, inputv, dim, smoo)
%inputv represent how many variable we have. Usually we 
%have two inputs after inputv: dim and smoo. If we have both inputs,
%inputv=[1,1];
% ; smoo strictly for display purposes
% ; motion.pro			4-17-00   Eric R. Weeks
% ;
% ; determines average displacement as a function of time for dt=1
% ; returns (t,dr,dx,dy)
% ;
% ; 6-13-05:  rewritten nearly from scratch, much better
% ;
% ; for more information see:
% ;       http://www.physics.emory.edu/~weeks/idl/motion.html
if inputv(1)==0; dim=2;end
if inputv(2)==0; smoo=10; end
mintime=1;  %in IDL it is automatically 1. We set it as 1 here.
ndat=length(tr(1,:));
dt=1;
totaltime=max(tr(:,ndat-1)); %in IDL it is max(tr(:,ndat-1),min=mintime). Here we omit min
result=double(zeros(totaltime+1,dim+2));
dx=getdx(tr,dt,1,dim);   %no problme for above
[Y, s]=sort(tr(:,ndat-1));  %matlab and IDL gives different result for s
[B,u]=unique(tr(s,ndat-1));
nu=length(u);
u=[0;u]; %in IDL it is -1 instead of 0
for i=2:nu+1
    tr0=tr(s(u(i-1)+1:u(i)),:);
    t0=tr0(1,ndat-1);
    dx0=dx(s(u(i-1)+1:u(i)),:);
    result(t0-mintime+1,1)=t0;  %set mintime=1
    w=find(dx0(:,dim+1)>-0.5);
    nw=length(w);
	if nw>0
		dx00=dx0(w,:);
		result(t0-mintime+1,3:dim+2)=sum(dx00(:,1:dim),1)/nw;
		result(t0-mintime+1,2) = nw;
    end
end

w=find(result(:,2)>0);
nw=length(w);
foo=result(w,:);
foo(:,2)=sqrt(sum(foo(:,3:dim+2).^2,2));
result(w,:)=foo; %we have the right result except that we haven't used rm_motion yet.
% temp=rm_motion(tr, result(w,:),inputv,dim, /just, smoo);
result=result(w,:);
