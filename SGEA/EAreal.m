function QWE=EAreal(Pop, EP, D, lu)
    % reproduction operation to produce children
    
    % mating selection is needed for most operators
    if rand <0.5
        ind1=EP(randi([1, size(EP,1)]),:); % choose one from Archive
        while true
            ind2_id= tournament_selection(2,1, Pop(:,end-1));% choose one from current population
            ind2=Pop(ind2_id,1:end-1);
            
            if ~all(ind1==ind2), break;end % make sure two selected inds are different
        end
        mating=[ind1; ind2];
    else
        while true
            mating_id=tournament_selection(2,2, Pop(:,end-1));
            mating =Pop(mating_id,1:end-1);

            if ~all(mating(1,:)==mating(2,:)), break;end % make sure two selected inds are different
        end
    end
    

    QWE = EAreal_sub(mating,D,lu); % this generates two children
    QWE = QWE(randi([1,2]),:);
end


function index = tournament_selection(K,N,F)
%TournamentSelection - Tournament selection
%
%   P = TournamentSelection(K,N,fitness) returns the indices
%   of N solutions by K-tournament selection based on their fitness values.
%   In each selection, the candidate having the minimum fitness1 value will
%   be selected.
%
%   Example:
%       P = TournamentSelection(2,100,FrontNo)
    [~,rank] = sortrows(F);
    [~,rank] = sort(rank);
    Parents  = randi(length(F),K,N);
    [~,best] = min(rank(Parents),[],1);
    index    = Parents(best+(0:N-1)*K);
end


function OffspringDec = EAreal_sub(ParentDec,D,lu)

%%%Parameters
proC=0.9;
disC=10;
proM=1;
disM=20;

[N, ] = size(ParentDec);

Parent1Dec = ParentDec(1:N/2,1:D);
Parent2Dec = ParentDec(N/2+1:end,1:D);
 %% Simulated binary crossover
beta = zeros(N/2,D);
mu   = rand(N/2,D);
beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
beta = beta.*(-1).^randi([0,1],N/2,D);
beta(rand(N/2,D)<0.5) = 1;
beta(repmat(rand(N/2,1)>proC,1,D)) = 1;
OffspringDec = [(Parent1Dec+Parent2Dec)/2+beta.*(Parent1Dec-Parent2Dec)/2
    (Parent1Dec+Parent2Dec)/2-beta.*(Parent1Dec-Parent2Dec)/2];


 %% Polynomial mutation
    Lower = repmat(lu(1,:),N,1);
    Upper = repmat(lu(2,:),N,1);
    Site  = rand(N,D) < proM/D;
    mu    = rand(N,D);
    temp  = Site & mu<=0.5;
    OffspringDec(temp) = OffspringDec(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                         (1-(OffspringDec(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    OffspringDec(temp) = OffspringDec(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                         (1-(Upper(temp)-OffspringDec(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
end