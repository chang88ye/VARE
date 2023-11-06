function ind=gaussian_mutate(ind,prob,lu)

%x=ind(D+1:Dim);

parDim=length(ind);

sigma=(lu(2,:)-lu(1,:))./20;

Q=normrnd(ind,sigma);

newparam=min(max(Q,lu(1,:)),lu(2,:));

for i=1:parDim
    if rand<prob
        ind(i)=newparam(i);
    else
        ind(i)=ind(i);
    end
end