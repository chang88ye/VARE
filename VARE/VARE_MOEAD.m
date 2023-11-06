function [Archives,Time, IGDs, HVs, SPs]=VARE_MOEAD(probID,N,runMax)
% <algorithm> 
% Vector Autoregressive Evolution for Dynamic Multi-Objective Optimization
%
% VARE with MOEA/D+DE as optimization engine

% Arguments:
%   INPUT:  probID  - specify the ID of the problem to be solved
%           N:      - specify population size
%           runMax  - specify maximum number of runs
%
%   OUTPUT: Archives - archived pop series for all runs
%           Time     - a vector of elapsed time for each run
%           IGDs     - a matrix of IGD values (runMax by the number of changes Ts)
%           HVs      - a matrix of HV values (runMax by the number of changes Ts)
%           SPs      - a a matrix of spacing values (runMax by the number of changes Ts)
% 
% Implemented by Dr. Shouyong Jiang on 4/1/2023.

% turn off warnings raised by PCA processes
warning('off', 'stats:pca:ColRankDefX');

% Parameter setting
[delta,nr] = deal(0.9,2);

Lag=5;   % time windows for prediction 

prob=problem(probID);
[varDim, objDim, bounds, func]=deal(prob.varDim, prob.objDim, prob.bounds, prob.func); % number of variables/objectives, variable bounds, and function handle of the problem

% get change severity/frequency and the past number of generations at which first change occurs

% method 1: get default values from problem definition
%[T0, taut, nt] =deal(instance.dynPar(1),instance.dynPar(2),instance.dynPar(3));

% method 2: get global values initialised in the main function
[nt, taut, T0]=dynamic_param(); 

Tg=T0+10*nt*taut; %set the termination genenerations: {T0 +10*nt*taut}

%% Generate the reference directions (general approach)
[W,N] = UniformPoint(N,objDim);
Tn = ceil(N/10); % neighbourhood size

%% Detect the neighbours of each solution
B = pdist2(W,W);
[~,B] = sort(B,2);
B = B(:,1:Tn);

%%
for run=1:runMax
    L=Lag;
    Lp=zeros(L,N);  % Times of prediction used by each search direction
    Ls=zeros(L,N);  % Times of success after applying prediction
    Lm=zeros(L,N);  % Times of success after applying hypermutation
    
    EP=[]; % external archive to save reference-wise pop before each change
    Archive=[];% archive to save pop before each change

    g=1; % generation counter
    %%%%%%%%%%%%%%%%%% Initialise population
    Qop=init_Pop(varDim,N,bounds);
    Pi=0.5*ones(1,N); % Probabiity that each search direction uses prediction
    
    %%%%%%% Evaluate initial population
    Pop= func(Qop, 0); % t=0
    Z = min(Pop(:, end-objDim+1:end),[],1);
    
    changeCounter=0; % to count the number of changes so far
    
    %%%%%%%%%%%%%%%%%%%%% Optimization loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tStart=tic;
    while g<=Tg 
        t=prob.update_time(g, taut, nt, T0); % update time instant
        checkvalue=checkChange(func, Pop,objDim, t, 0.1); % 10% population for change detection
        disp(['Evolve at generation ',num2str(g), ', time t= ', num2str(t), '; ----->']);

        if checkvalue==1  %% when a change is detected
            changeCounter=changeCounter+1;
            t_prev=prob.update_time(g-1, taut, nt, T0);
            
            [~, Ei]=Associate(Pop(:,end-objDim+1:end),W);
            EP=[EP,Pop(Ei,:)]; % create time-series for each w
            
           [IGD,SP,HV]=Performance(Pop(:, end-objDim+1:end),prob,t_prev);% calculate performance metrics
            
            IGDt(changeCounter)=IGD;
            SPt(changeCounter)=SP;
            HVt(changeCounter)=HV;

            Archive=[Archive;Pop];
    
            [Pop, Offspring, Predicted]=ChangeResponse(EP, Pop,Pi,Lag, W, prob, t);%respond to the change
            
            % Update the ideal point
            Z = min([Z; Offspring(:,varDim+1:end)], [], 1);

            
            % Update the solutions in P by Tchebycheff approach
            g_old = max(abs(Pop(:,end-objDim+1:end)-Z).*W,[],2);
            g_new = max(abs(Offspring(:,end-objDim+1:end)-Z).*W,[],2);
            Chosen =find(g_old>g_new);
            Pop(Chosen,:)=Offspring(Chosen,:);

            idx=mod(changeCounter-1, L)+1;
            Lp(idx, :)=0;
            Lp(idx, Predicted)=1;
            
            % update pi which is the probability to choose prediction as change response strategy
            pred_succeed=intersect(Chosen,Predicted); % predicted individuals that are sucessfully kept in the new pop
            mut_succeed =setdiff(Chosen,pred_succeed); % mutated individuals that are sucessfully kept in the new pop
            
            Ls(idx, :)=0;
            Ls(idx, pred_succeed)=1;
            Lm(idx, :)=0;
            Lm(idx, mut_succeed)=1;
            
            % update only when there are enough lags 
            if changeCounter>=Lag      
                  Pp=sum(Ls,1)./(1e-6+sum(Lp,1));
                  Pm=sum(Lm,1)./(1e-6+L-sum(Lp,1));
                  Pi=Pp./(Pp+Pm);
                  Pi(isnan(Pi))=0.5;
                  Pi= min(max(0.1, Pi), 0.9);
            end
            
        else
        
            for i=1:N

                % Choose the parents
                if rand < delta
                    P = B(i,randperm(size(B,2)));
                else
                    P = randperm(N);
                end

                % Generate an offspring by DE
                QWE = DE_operator(bounds,Pop([i, P(1:2)],1:varDim));
                Qop = checkbound(QWE,bounds);
                Offspring = func(Qop, t); %% evaluate the new population

                % Update the ideal point
                Z = min([Z; Offspring(:,varDim+1:end)], [], 1);

                % Update the solutions in P by Tchebycheff approach
                g_old = max(abs(Pop(P,varDim+1:end)-repmat(Z,length(P),1)).*W(P,:),[],2);
                g_new = max(repmat(abs(Offspring(varDim+1:end)-Z),length(P),1).*W(P,:),[],2);

                idx=P(find(g_old>=g_new,nr));
                Pop(idx,:) = repmat(Offspring, length(idx),1);
            end
        end
        
        
        g=g+1;
    end
    %%%%%%
    tEnd=toc(tStart);
    
    Time(run)=tEnd;
    IGDs(run,:)=IGDt;
    SPs(run,:)=SPt;
    HVs(run,:)=HVt;
    Archives{run}=Archive;
end
end

function [P_new, OffSpring, to_predict]=ChangeResponse(EP, Pop, Prob, Lag, W, instance,t)
N=size(W, 1);
dim=size(Pop, 2);

[D, lu, func]=deal(instance.varDim, instance.bounds, instance.func);

P_cur=EP(:,end-dim+1:end);%% curr population which has done association with W

% re-evaluation
P_new= func(P_cur(:, 1:D), t); 

% identification of individuals for each reference direction
[~, Ei_new]=Associate(P_new(:,D+1:end),W);

delta_F=mean(abs((P_new(Ei_new,D+1:end)-P_cur(:,D+1:end))./P_cur(:,D+1:end)), 'all'); % delta_F: amount of change in objective space

delta_X=mean(abs((P_new(Ei_new,1:D)-P_cur(:,1:D))./P_cur(:,1:D)), 'all'); % delta_X: amount of change in decision space

OffSpring=P_cur(Ei_new,1:D); % population - values of decision variables only

% apply hypermutation for ealy changes
HM=rand(1,N)>Prob; % individuals that perform hypermutation
to_predict=find(~HM); % individuals to undergo prediction
       
predict_OK=true;
if size(EP,2)>=dim*4*Lag % varm requires at least p^2*lag data samples (we assume there are p>=2 variables after PCA for any problem)
    for i = to_predict
        input_data = reshape (EP(i,:), dim, [])';
        input_data = input_data(:,1:D);
        % PCA processing: https://stats.stackexchange.com/questions/229092/how-to-reverse-pca-and-reconstruct-original-variables-from-several-principal-com
        [coeff,score,~,~,explained,mu] = pca(input_data);
        nComp = find(cumsum(explained) >= 80,1); % find the first k largest PCs that explain the data >80%
        Y     = score(:,1:nComp);

        try
			mdl     = bayesvarm(nComp,Lag);
            EstMdl  = estimate(mdl,Y, 'Display', 0);
            x = forecast(EstMdl,1,Y);
            OffSpring(i,:)=x(:,1:nComp) * coeff(:,1:nComp)'+mu;
            if ~isreal(OffSpring(i,:))
                predict_OK=false;
                break;
            end
        catch % most likely not enough data to build varm models
            predict_OK=false;
            break
        end
    end
else
    predict_OK=false;
end

if ~predict_OK
    if size(EP,2)>dim % data has no less than 2 lags and less than 4*L lags, perform one-step linear prediction
        for i = to_predict
            input_data   = reshape (EP(i,:), dim, [])';
            input_data = input_data(:,1:D);
            [x_cur, x_pre]=deal(input_data(end,:), input_data(end-1,:));
            OffSpring(i,:)=x_cur+(x_cur - x_pre)+ randn(1,1)*norm(x_cur - x_pre) / (2*D^0.5);
        end
    else
        HM = true(1, N);
        to_predict=[];
    end
end

%%%%% the following statistic test is deactivated as PS change often occurs in dynamic environments
% t-test to check whether there is significant change is PS
% tempDiff=abs(P_new(Ei_new,1:D)-P_cur(:,1:D)); % a N by D matrix
% h = ttest(tempDiff(:),0,"Tail","right"); %test the null hypothesis that the data comes from a population with mean equal to 0, against the alternative that the mean is greater than 0.
% if h >0
%     disM=20*max(exp(-(delta_F+delta_X)),0.1); %mean([exp(-delta_F), exp(-delta_X)]);
% else
%     disM=20;
% end

% do hypermutation
disM=20*max(exp(-(delta_F+delta_X)),0.1); %mean([exp(-delta_F), exp(-delta_X)]);

OffSpring(HM,:)=AdaptiveMutation(OffSpring(HM,:), disM, lu, sum(HM), D);

Qop=checkbound(OffSpring,lu);

% evaluate population in new environment
OffSpring=func(Qop,t);
end

function [dist, Ei]=Associate(PopObj,W)
QObj    = ObjectiveNormalization(PopObj);
cosine = 1-pdist2(QObj,W, 'cosine');
[cosVal, Ei] = max(cosine, [], 1); % find the closest ind for each w_i
dist= vecnorm(QObj(Ei,:),2,2).*sin(acos(cosVal))'; %1-cosVal; %
end

function PopDec=AdaptiveMutation(PopDec, disM, lu, N, D)
%% Polynomial mutation
Lower = repmat(lu(1,:),N,1);
Upper = repmat(lu(2,:),N,1);
Site  = rand(N,D) < 1/D;
mu    = rand(N,D);
temp  = Site & mu<=0.5;
PopDec(temp) = PopDec(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
    (1-(PopDec(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
temp = Site & mu>0.5;
PopDec(temp) = PopDec(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
    (1-(Upper(temp)-PopDec(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
end
 

function OffspringDec = DE_operator(bounds,ParentDecs)
    [CR,F,proM,disM] = deal(1,0.5,1,20);
    ParentDecs    = ParentDecs([1:end,1:ceil(end/3)*3-end],:);
    [N,D]     = size(ParentDecs);

    %% Differental evolution
    Parent1Dec   = ParentDecs(1:N/3,:);
    Parent2Dec   = ParentDecs(N/3+1:N/3*2,:);
    Parent3Dec   = ParentDecs(N/3*2+1:end,:);
    OffspringDec = Parent1Dec;
    Site = rand(N/3,D) < CR;
    OffspringDec(Site) = OffspringDec(Site) + F*(Parent2Dec(Site)-Parent3Dec(Site));

    %% Polynomial mutation
    Lower = repmat(bounds(1,:),N/3,1);
    Upper = repmat(bounds(2,:),N/3,1);
    Site  = rand(N/3,D) < proM/D;
    mu    = rand(N/3,D);
    temp  = Site & mu<=0.5;
    OffspringDec(temp) = OffspringDec(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                         (1-(OffspringDec(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    OffspringDec(temp) = OffspringDec(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                         (1-(Upper(temp)-OffspringDec(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
end