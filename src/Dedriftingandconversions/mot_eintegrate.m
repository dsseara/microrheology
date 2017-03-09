function result=mot_eintegrate(data)
%written by Yongxiang Gao on Nov28 2005
% result=reform(data)
result=data;
% result(1)=data(1);
% result(1)=0.0;
n=length(data);
for i=2:n
    result(i)=result(i-1)+data(i);
end
% result(1)=data(1);