function [subp,neighbour_indexMatrix]=init_weights(popsize,niche,ObjNum)
%输出 权重矩阵 及 每个个体对应的邻域个体指标
subp=[];%%% 权重矩阵
neighbour_indexMatrix=[]; %%每个个体对应的邻域个体指标
for i=1:popsize
    if ObjNum==2
        Weight=linspace(i,popsize-i,ObjNum);
        weight=Weight./popsize;
%         weight=zeros(1,2);
%         weight(1)=i/popsize;
%         weight(2)=(popsize-i)/popsize;
        subp=[subp;weight];
    elseif ObjNum==3
        Weight=linspace(i,popsize-i,ObjNum);
        weight=Weight./popsize;
        subp=[subp;weight];
    end
end

leng=length(subp);%矩阵中较大的指标
distanceMatrix=zeros(leng,leng);

for i=1:leng
    for j=i+1:leng
        A=subp(i,:);
        B=subp(j,:);
        distanceMatrix(i,j)=(A-B)*(A-B)'; %求与任意个体之间的距离
        distanceMatrix(j,i)= distanceMatrix(i,j);%作对称
    end
    [s,index]=sort(distanceMatrix(i,:));
    neighbour_indexMatrix(i,:)=index(1:niche);
end
        
         
        
        
        
        
        
        