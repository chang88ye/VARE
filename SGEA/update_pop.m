function [pop, EP]=update_pop(ind, pop, EP, ObjNum)
% update current population and external population that act as archive
% note that pop has the struture [dec_vector, obj_vector, scalar_fitness]
N=size(pop,1);
indObj=ind(1,end-ObjNum+1:end);
popObj=pop(:,end-ObjNum:end-1);
dominated=0;
noncomp=0;
for j = 1 : N
    k = any(indObj<popObj(j,:)) - any(indObj>popObj(j,:));
    if k == 1
        pop(j,end)=pop(j,end)+1;
    elseif k == -1
        dominated = dominated +1;
    else
        noncomp = noncomp +1;
    end
end
% identify the worst member in pop. If it is no better than ind, replace it with ind
[value, worst]=max(pop(:,end));
if value>=dominated
   if noncomp==N
       worst = randi([1,N]);
   end
   pop(worst,:) = [ind, dominated];
end

if dominated <1 % ind is nondominated by pop
    EPObj=EP(:,end-ObjNum+1:end);
    dominated=false(1,size(EP,1));
    for j = 1 : size(EP,1)
        k = any(indObj<EPObj(j,:)) - any(indObj>EPObj(j,:));
        if k == 1
            dominated(j)=true;
        elseif k == -1
            break;
        end
    end
    EP(dominated,:)=[];
    EP=[EP; ind];
end
end