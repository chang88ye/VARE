function Score = IGD(PopObj,PF)
% <metric> <min>
% Inverted generational distance

if nargin<2
    Score=[];
else
    Distance = min(pdist2(PF,PopObj),[],2);
    Score    = mean(Distance);
end

end