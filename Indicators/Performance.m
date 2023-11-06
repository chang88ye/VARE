% Performance function
function [IGDScore,SPScore,HVScore]=Performance(PopObj,prob,t)

    % get the true PF and approximation
    Pareto=generatePF(prob,t);
    
    % normalisation
    if ~isempty(Pareto)
        Fmax   = max(Pareto, [], 1);
        Fmin   = min(Pareto, [], 1);
    
        % degenerate case where Fmax=Fmin, for which Fmax is used.
        Fmax_min=Fmax - Fmin;
        Fmax_min(Fmax_min<1e-3)=max(Fmax(Fmax_min<1e-3), 1e-3); % just in case Fmax=0
    
        Pareto = (Pareto - Fmin)./Fmax_min;
        PopObj = (PopObj - Fmin)./Fmax_min;

        %%%%%%%HV score
        HVScore = HV(PopObj,Pareto);
        %%%%%%%IGD score
        IGDScore=IGD(PopObj,Pareto);
        %%%%%%%SP score
        SPScore=Spacing(PopObj,Pareto);
    else
        HVScore = HV(PopObj);
        %%%%%%%IGD score
        IGDScore=-1;
        %%%%%%%SP score
        SPScore=Spacing(PopObj);
    end

end

