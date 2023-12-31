function Population = EnvironmentalSelection(Population,N,D,ObjNum)
% The environmental selection of RM-MEDA

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population(:,D+1:D+ObjNum),N);    
    Next = FrontNo < MaxFNo; %%0,1

    %% Delete the solutions in the last front by crowding distance
    Last = find(FrontNo==MaxFNo);

    while length(Last) > N - sum(Next)
        [~,worst]   = min(CrowdingDistance(Population(Last,D+1:D+ObjNum)));
        Last(worst) = [];
    end
    Next(Last) = true;

    % Population for next generation
    
    Population = Population(Next,1:D+ObjNum); 

end

function CrowdDis = CrowdingDistance(PopObj)
% Calculate the crowding distance of each solution in the same front

    [N,M]    = size(PopObj);
    
    CrowdDis = zeros(1,N);
    Fmax     = max(PopObj,[],1);
    Fmin     = min(PopObj,[],1);
    for i = 1 : M
        [~,rank] = sortrows(PopObj(:,i));
        CrowdDis(rank(1))   = inf;
        CrowdDis(rank(end)) = inf;
        for j = 2 : N-1
            CrowdDis(rank(j)) = CrowdDis(rank(j))+(PopObj(rank(j+1),i)-PopObj(rank(j-1),i))/(Fmax(i)-Fmin(i));
        end
    end
end