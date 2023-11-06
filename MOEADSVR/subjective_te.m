function obj=subjective_te(weight,ind,idealpoint,D)
%输出的是一个列向量或者一个数值。
%weight是20*2的
%只适合于一个个体
s=size(weight,1);  %niche
indsize=size(ind,1); %可能是niche个，也可能是一个

weight((weight==0))=0.00001;

if indsize==s 
   % size(ind(:,D+1:Dim))
  %  idealpoint
     %size(idealpoint(ones(1,indsize),:))
%ind(:,D+1,Dim)-idealpoint(ones(1,indsize),:)
    part2=abs(ind(:,D+1:end)-idealpoint(ones(1,indsize),:));
    obj=max(weight.*part2,[],2);%%要按行求最大，输出一个列向量
elseif indsize==1   %当只有一个个体地情况下
    part2=abs(ind(D+1:end)-idealpoint);
    
    obj=max(weight.*part2(ones(1,s),:),[],2);
else
    error('反正是不对');
end