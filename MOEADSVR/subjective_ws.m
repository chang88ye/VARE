function obj=subjective_ws(weight,ind,D)
 s=size(weight,1);
indsize=size(ind,1); 
weight((weight==0))=0.00001;

if indsize==1
    part2=ind(D+1:end);
    obj=sum(weight.*part2(ones(1,s),:),2);
else
    obj=sum(weight.*ind(:,D+1:end),2);
end
