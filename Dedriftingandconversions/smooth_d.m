function r=smooth_d(oridata,w)
%simplely smooth data linearly. Corresponding to the smooth function in
%IDL. Written by Yongxiang Gao Nov28, 2005.
n=length(oridata);
for i=1:n
    if i>=(w+1)/2 && i<=n-(w-1)/2
        if mod(w,2)==1
        r(i)=sum(oridata(i+1-round((w+1)/2):i+w-fix((w+1)/2)))/w;
        else
        r(i)=sum(oridata(i+1-round((w+1)/2):i+w-fix(w/2)))/(w+1);  
        end
    else
        r(i)=oridata(i);
    end
end