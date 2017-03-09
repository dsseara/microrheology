function rdf=dedrift(lub, dr)
%written by Yongxiang Gao on June 2, 2006
%to remove drift
for i=1:length(dr(:,1));
     id=find(lub(:,3)==i+1);
     lub(id,1)=lub(id,1)-dr(i,1);
     lub(id,2)=lub(id,2)-dr(i,2);
end
    rdf=lub;