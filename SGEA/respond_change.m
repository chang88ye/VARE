function [Pop, EP, C0]=respond_change(ins, EP, Pop, C0, eta, t)
% change response consisting of 50% current population, eta% random
% immigrants, and predicted individuals

    N=size(Pop,1);
    D=ins.varDim;
    
    PopObjs=Pop(:,D+1:D+ins.objDim);
    
    C1=mean(EP(:,1:D),1);% centroid of archive
    
    %% method 1
    %
    % select extreme solutions to guid FF selection method
    [keep,selected]=filter_boundary(PopObjs);
    
    % calculate the nearest members in selected for each of the unselected
    dist=pdist2(PopObjs(selected,:),PopObjs(keep,:),'euclidean','smallest',1);
    
    % select half of population using the-farthest-the-first method
    idx_undel=find(keep);
    for j=1:floor(N/2)-length(selected)
        [~, ind]=max(dist);
        selected=[selected, idx_undel(ind)];
        dist=min(dist,pdist2(PopObjs(idx_undel(ind),:),PopObjs(keep,:),'euclidean'));   
    end
    %%
    % method 2
    %{
    [~,ind1]=min(Pop(:,D+1));
    [~,ind2]=max(Pop(:,D+1));
    selected=unique([ind1,ind2]);
    
    
    toSel=setdiff(1:N, selected);
    
    % calculate the nearest members in selected for each of the unselected
    dist=pdist2(PopObjs(selected,:),PopObjs(toSel,:),'euclidean','smallest',1);
    
    % select half of population using the-farthest-the-first method
    for j=1:floor(N/2)
        [~, ind]=max(dist);
        selected=[selected, toSel(ind)];
        dist=min(dist,pdist2(PopObjs(toSel(ind),:),PopObjs(toSel,:),'euclidean'));   
    end
    %}
    
    % split population
    remainded=setdiff(1:N, selected);
    immigrants=rand(1,length(remainded))<2*eta;
    remainded=remainded(~immigrants);
    
    % create random immigrants
    Qop_imm=init_population(ins,sum(immigrants));
    Qop_sel=Pop(selected,1:D);
    
    Qop1=repair([Qop_sel;Qop_imm],ins.bounds);
    Qop1=ins.func(Qop1(:, 1:D),t);
    
    % identify nondominated inds
    [~, nondomIn] = calStrength(Qop1(:, D+1:end));
    
    Cn=mean(Qop1(nondomIn,1:D),1); % nondominated centroid of responded pop at new envionment
    C2=mean(Qop1(:,1:D),1); % centroid of responded pop
    
    %prediction
    if ~isempty(C0)
        % determine stepsize (C0 - prev archive centroid, C1 -current archive centroid) 
        stepsize = norm(C0-C1);
    
        Qop2 = Pop(remainded,1:D)+stepsize*normalize(Cn-C2,'norm') + ...
            randn(1, D)*stepsize/sqrt(4*D);
    else % random immigrants if it is the first change
        Qop2 =init_population(ins,length(remainded));
    end
    
    Qop2=repair(Qop2,ins.bounds);
    Qop2=ins.func(Qop2(:,1:D),t);
    
    Pop=[Qop1; Qop2];
    C0=C1;
    
    % assign fitness for pop and update EP
    [Pop,EP]=environmental_selection([], Pop, ins.objDim);
end

