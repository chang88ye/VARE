function [Pop, Fevals]=Tr_IPG_response(TotalArchive,lu,ObjNum,popsize,changetime,func,t)

%% Note: this kind of repsonse in each environmental change has three rounds of fitness evaluation
% The total FEs in each response is 3*sampleN =3*popsize (if sampleN=popsize), 
% much higher than other type of response mechanisms thathave popsize FEs.
%% 

sampleN=popsize;

Fevals=0;

[~,Dim]=size(TotalArchive);
D=Dim - ObjNum;

pop_current=TotalArchive((changetime-1)*popsize+1:changetime*popsize,:);%%当前代种群

tempParticle=init_Pop(D,sampleN,lu);
QSFtemp=func(tempParticle,t-1); % first round of evaluation: sampleN FEs
Fs=QSFtemp(:,D+1:Dim);   %%%这块只需要取相应的目标函数值就可以了
Fevals = Fevals+sampleN;

ttempParticle=init_Pop(D,sampleN,lu);
QSFttemp=func(ttempParticle,t); % second random of evaluation: sampleN FEs
Fa=QSFttemp(D+1:Dim);%%%只是取相应的函数值
Fevals=Fevals+sampleN;


% Find the latent space of domain adaptation
mu = 0.5;
lambda = 'unused';
dim = 20;           % Deduced dimension
kind = 'Gaussian';  % The dimension of Gaussian Kernel feature space is inifinite, so the deduced dimension can be 20.
p1 = 1;
p2 = 'unused';
p3 = 'unused';
W = getW(Fs', Fa', mu, lambda, dim, kind, p1, p2, p3);  %%%TCA计算W  %%%到这都没有问题

POF_deduced = getNewY(Fs', Fa', pop_current(:,D+1:Dim)', W, kind, p1, p2, p3);

% dis_px = @(p, x)sum((getNewY(Fs', Fa', obj_func(x,D,ObjNum,Dim,problem,t)', W, kind, p1, p2, p3) - p).^2);
dis_px = @(p, x)sum((getNewY(Fs', Fa', obj_func(x,D, func,t)', W, kind, p1, p2, p3) - p).^2);
initn = size(POF_deduced, 2);

for i = 1:initn % search in the decision space
    start=lu(1,:)+(lu(2,:)-lu(1,:)).*rand(1,D);
%     indiv= fmincon(@(x)dis_px(POF_deduced(:,i), x), start, ...
%         [], [], [], [], lu(1,:), lu(2,:), 'Display', 'off');
    indiv= fmincon(@(x)dis_px(POF_deduced(:,i), x), start, ...
        [], [], [], [], lu(1,:), lu(2,:), [], optimset('display', 'off'));


    % evaluate new individual
    Pop(i,:)=func(indiv,t); % third round of evaluation
    Fevals=Fevals+1;
end

end

% helper function
function objs=obj_func(x,D,func,t)
xy=func(x,t);
objs=xy(:,D+1:end);
end

% function objs=obj_func(x,D,ObjNum,Dim,problem,t)
% xy=benchmark_func(x,D,ObjNum,Dim,problem,t);
% objs=xy(:,D+1:end);
% end
%             
            
            