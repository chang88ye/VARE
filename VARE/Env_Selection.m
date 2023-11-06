function [Population, Choose] = Env_Selection(stop,Population,Ei,FV,N)
% The environmental selection of SPEA/R

Choose = [];
while length(Choose) < N
    H = [];
    for i = unique(Ei)
        if i > 0
            Local = find(Ei==i);
            [~,q] = min(FV(Local));
            H     = [H,Local(q)];
        end
    end
    if  stop && any(FV(H)<1)
        H=H(FV(H)<1);
    end
    Ei(H)=-1;
    if length(Choose) + length(H) < N
        Choose = [Choose,H];
    else
        [~,rank] = sort(FV(H));
        Choose   = [Choose,H(rank(1:N-length(Choose)))];
    end
end
Population = Population(Choose,:);

end