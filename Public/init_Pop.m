function Pop=init_Pop(D,popsize,lu)
    % random solutions in the range of [lb, ub]
    randarray=rand(popsize,D);
    
    lowend=lu(1,:);
    upend=lu(2,:);
    span=upend-lowend;
    
    Pop=randarray.*(span(ones(1,popsize),1:D))+lowend(ones(1,popsize),1:D);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
