function [Archives,Time, IGDs, HVs, SPs]=VARE(probID,N,runMax)
% <algorithm> 
% Vector Autoregressive Evolution for Dynamic Multi-Objective Optimization

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
Lag=5;   % time windows for prediction 

prob=problem(probID);
disp(['==================== Solving problem ',prob.name, ' ======================']);
[varDim, objDim, bounds, func]=deal(prob.varDim, prob.objDim, prob.bounds, prob.func); % number of variables/objectives, variable bounds, and function handle of the problem

% get change severity/frequency and the past number of generations at which first change occurs

% method 1: get default values from problem definition
%[T0, taut, nt] =deal(instance.dynPar(1),instance.dynPar(2),instance.dynPar(3));

% method 2: get global values initialised in the main function
[nt, taut, T0]=dynamic_param(); 

Tg=T0+10*nt*taut; %set the termination genenerations: {T0 +10*nt*taut}

%% Generate the reference directions (general approach)
[W,N] = UniformPoint(N,objDim);
% Largest acute angle between two neighbouring reference directions
cosine = 1 - pdist2(W,W,'cosine');
cosine(logical(eye(length(cosine)))) = 0;
theta  = max(min(acos(cosine),[],2));

%%
for run=1:runMax
    EP=[]; % external archive to save reference-wise pop before each change
    Archive=[];% archive to save pop before each change

    g=1; % generation counter
    %%%%%%%%%%%%%%%%%% Initialise population
    Qop=init_Pop(varDim,N,bounds);
    Pi=0.5*ones(1,N); % Probabiity that each search direction uses prediction

    L=Lag;
    Lp=zeros(L,N);  % Times of prediction used by each search direction
    Ls=zeros(L,N);  % Times of success after applying prediction
    Lm=zeros(L,N);  % Times of success after applying hypermutation
    
    %%%%%%% Evaluate initial population
    Pop= func(Qop, 0); % t=0
    
    changeCounter=0; % to count the number of changes so far
    
    %%%%%%%%%%%%%%%%%%%%% Optimization loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tStart=tic;
    while g<=Tg 
        t=prob.update_time(g, taut, nt, T0); % update time instant
        checkvalue=checkChange(func, Pop,objDim, t, 0.1); % 10% population for change detection
        disp(['Evolve at generation ',num2str(g),  ' for ', prob.name, ', time t= ', num2str(t), '; ----->']);

        if checkvalue==1  %% when a change is detected
            changeCounter=changeCounter+1;
            t_prev=prob.update_time(g-1, taut, nt, T0);
            
            [~, Ei]=Associate(Pop(:,end-objDim+1:end),W);
            EP=[EP,Pop(Ei,:)]; % create time-series for each w
            
            [IGD,SP,HV]= Performance(Pop(:, end-objDim+1:end),prob,t_prev);%%% calculate performance metrics; 
            
            IGDt(changeCounter)=IGD;
            SPt(changeCounter)=SP;
            HVt(changeCounter)=HV;

            Archive=[Archive;Pop];
    
            [Pop, Offspring, Predicted]=ChangeResponse(EP, Pop,Pi,Lag, W, prob, t);%respond to the change

            idx=mod(changeCounter-1, L)+1;
            Lp(idx, :)=0;
            Lp(idx, Predicted)=1;
        else
            % Staticc Evolutionary Multiobjective Optimization Process
            QWE = RMMEAD_operator(Pop,varDim,objDim);%%%%% generate new population suing RM-MEDA operator
            Qop=checkbound(QWE,bounds);
            Offspring= func(Qop, t); %% evaluate the new population
        end
        
        % environmental selection using diversity-based sorting
        QObj       = ObjectiveNormalization([Pop(:,varDim+1:end); Offspring(:,varDim+1:end)]);
        Dist       = pdist2(W,QObj,'cosine'); % pair-wise cosine distance

        [Angle, Ei]= min(acos(1-Dist), [], 1);
        FV         = FitnessAssignment(Ei,QObj,Angle,theta);
        [Pop, Chosen] = Env_Selection(g==Tg,[Pop; Offspring],Ei,FV,N);
        
        % update pi which is the probability to choose prediction as change response strategy
        if checkvalue==1
            Chosen=Chosen(Chosen>N)-N; % indices of selected offspring
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
            EstMdl  = estimate(mdl,Y,'Display',0);
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
 