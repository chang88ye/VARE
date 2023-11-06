function popDecs=init_population(problem,popsize)
% Uniform sampling of decision space

% Input:
%   problem -problem struct
%   popsize -population size

% Output:
%   popDecs - population in decision space

    lower=problem.bounds(1,:);
    upper=problem.bounds(2,:);
    
    randarray=rand(popsize,problem.varDim);
    
    popDecs=randarray.*(upper-lower)+lower;
end