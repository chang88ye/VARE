function [f,p]=NondomSort(u,D,ObjNum,Dim)
% [f,p]=NondomSort(u,D,ObjNum,Dim)
% nondominated sorting to identify nondominated individuals
% p - indices of nondominated individuals
% f - nondominated individuals
m=size(u,1);
p=[];    
for i=1:m
    num=0;
    for j=1:m
        dom_less=0;
        dom_equal=0;
        dom_more=0;
        for k=1:ObjNum
            if u(i,D+k)<u(j,D+k)          
                dom_less=dom_less+1;             
            elseif u(i,D+k)==u(j,D+k)
                dom_equal=dom_equal+1;
            else
                dom_more=dom_more+1;
            end
        end
        if dom_less==0&&dom_equal~=ObjNum
            num=num+1;                        
        end
    end 
    if num==0                                  
       p=[p i];                             
    end
end
f=u(p,1:Dim);
