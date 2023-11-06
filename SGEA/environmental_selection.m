function [pop,EP]=environmental_selection(parents, offspring, ObjNum)
% pop: parents; offspring: children.
% population trucation by k-nearest neighbor based on sde-distance 
N=size(offspring,1);
comb_pop=unique([parents;offspring],'rows');
comObjs=comb_pop(:,end-ObjNum+1:end);
[strengths, nondominated]=calStrength(comObjs);

fitness=strengths(nondominated);
if length(strengths)>=N
    [~,selected]=mink(strengths, N);
    
    if length(nondominated)<N
        pop=[comb_pop(selected,:),strengths(selected)];
        EP=comb_pop(nondominated,:);
    elseif length(nondominated)==N
        pop=[comb_pop(nondominated,:), fitness];
        EP=comb_pop(nondominated,:);
    else
        tmp=comb_pop(nondominated,:);
        deleted=truncation(tmp(:,end-ObjNum+1:end), length(nondominated)-N);
        pop=[tmp(~deleted,:), fitness(~deleted)];
        EP=tmp(~deleted,:);
    end
else
    EP=comb_pop(nondominated,:);
    [fit,idx] = mink(strengths, N);

    com_size=size(comb_pop,1);
    if com_size<N % ensure popsize is maintained
        k=randperm(com_size,N-com_size);
        pop=[comb_pop, strengths; 
            comb_pop(k, :), strengths(k)];
    else
        pop=[comb_pop(idx,:), fit];
    end

end
if size(pop,1)<100
    a=1;
end
end

function Del = truncation(PopObj,K)
% Select part of the solutions by truncation

    N = size(PopObj,1);
    
    %% Calculate the shifted distance between each two solutions
    Distance = inf(N);
    for i = 1 : N
        SPopObj = max(PopObj,repmat(PopObj(i,:),N,1));
        for j = [1:i-1,i+1:N]
            Distance(i,j) = norm(PopObj(i,:)-SPopObj(j,:));
        end
    end
    
    %% Truncation
    Del = false(1,N);
    while sum(Del) < K
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end