function indPop=genetic_operate(neighbour_indexMatrix,Pop,index,F,CR,lu,D)
%%%%取出对应的邻域个体对应指标
neighbourindex=neighbour_indexMatrix(index,:);

nsize=length(neighbourindex);
%%%%%找三个不同的个体
r(1)=neighbourindex(ceil(rand*nsize));
while r(1)==index
    r(1)=neighbourindex(ceil(rand*nsize));
end

r(2)=neighbourindex(ceil(rand*nsize));
while r(2)==index||r(2)==r(1)
    r(2)=neighbourindex(ceil(rand*nsize));
end

r(3)=neighbourindex(ceil(rand*nsize));
while r(3)==index||r(3)==r(1)||r(3)==r(2)
    r(3)=neighbourindex(ceil(rand*nsize));
end

%找出这些指标对应的个体
Pop(r(1),:);
Pop(r(2),:);
Pop(r(3),:);

v=Pop(r(1),1:D)+F*(Pop(r(2),1:D)-Pop(r(3),1:D));

for j=1:D
    if rand<CR
        newpoint(j)=v(j);
    else
        newpoint(j)=Pop(index,j);
    end
end
%%边界处理

newpoint=max(newpoint,lu(1,:));
newpoint=min(newpoint,lu(2,:));
 

%%%个体基因变异操作
indPop=gaussian_mutate(newpoint,1/2,lu);






