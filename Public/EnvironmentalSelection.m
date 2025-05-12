function [Population, FrontNo, CrowdDis] = EnvironmentalSelection(Population,N,D,ObjNum)
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

    % %% Delete the solutions in the last front by crowding distance
    % Last = find(FrontNo==MaxFNo);
    % 
    % while length(Last) > N - sum(Next)
    %     [~,worst]   = min(CrowdingDistance(Population(Last,D+1:D+ObjNum)));
    %     Last(worst) = [];
    % end
    % Next(Last) = true;
    % 
    % % Population for next generation
    % 
    % Population = Population(Next,1:D+ObjNum); 

    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(Population(:,D+1:D+ObjNum),FrontNo);
    
    %% Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
    
    %% Population for next generation
    Population = Population(Next,:);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);

end

function CrowdDis = CrowdingDistance(PopObj,FrontNo)
% Calculate the crowding distance of each solution front by front

    [N,M]    = size(PopObj);
    CrowdDis = zeros(1,N);
    Fronts   = setdiff(unique(FrontNo),inf);
    for f = 1 : length(Fronts)
        Front = find(FrontNo==Fronts(f));
        Fmax  = max(PopObj(Front,:),[],1);
        Fmin  = min(PopObj(Front,:),[],1);
        for i = 1 : M
            [~,Rank] = sortrows(PopObj(Front,i));
            CrowdDis(Front(Rank(1)))   = inf;
            CrowdDis(Front(Rank(end))) = inf;
            for j = 2 : length(Front)-1
                CrowdDis(Front(Rank(j))) = CrowdDis(Front(Rank(j)))+(PopObj(Front(Rank(j+1)),i)-PopObj(Front(Rank(j-1)),i))/(Fmax(i)-Fmin(i));
            end
        end
    end
end