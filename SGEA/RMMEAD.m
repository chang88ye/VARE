function OffspringDec = RMMEAD(Model,probability, Num, D, M)

%% Reproduction
OffspringDec = zeros(Num,D);
% Generate new trial solutions one by one
for i = 1 : Num
    % Select one cluster by Roulette-wheel selection
    k = find(rand<=probability,1);
    % Generate one offspring
    if ~isempty(Model(k).eVector)
        lower = Model(k).a - 0.25*(Model(k).b-Model(k).a);
        upper = Model(k).b + 0.25*(Model(k).b-Model(k).a);
        trial = rand(1,M-1).*(upper-lower) + lower;
        sigma = sum(abs(Model(k).eValue(M:D)))/(D-M+1);
        OffspringDec(i,:) = Model(k).mean + trial*Model(k).eVector(:,1:M-1)' + randn(1,D)*sqrt(sigma);
    else
        OffspringDec(i,:) = Model(k).mean + randn(1,D);
    end
end

