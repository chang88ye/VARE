function [Archives,Time, IGDs, HVs, SPs]=MOEADSVR(probID,N,runmax)
% MOEA/D-SVR for Dynamic Multi-objective Optimization
%
% Modified by Dr. Shouyong Jiang on 4/1/2023

%%%%%%%%Algorithmic parameters%%%%%%%%%%%%%%%%%%%%
niche=floor(N/5);% size of neighbourhood
F=0.5;
CR=0.5;


prob=problem(probID);
[varDim, objDim, bounds, func]=deal(prob.varDim, prob.objDim, prob.bounds, prob.func); % number of variables/objectives, variable bounds, and function handle of the problem

% get change severity/frequency and the past number of generations at which first change occurs

% method 1: get default values from problem definition
%[T0, taut, nt] =deal(instance.dynPar(1),instance.dynPar(2),instance.dynPar(3));

% method 2: get global values initialised in the main function
[nt, taut, T0]=dynamic_param(); 

Tg=T0+10*nt*taut; %{T0 +10*nt*taut}

% %%---DHTSP setting start --- 
% %%comment out this block if other problems are tested
% 
% % Tmax=48; % 48 hours optimisation period
% Tg=960; % total number of generations
% M=Tg/taut; % number of changes 32 (taut=30) 48 (taut=20), 60 (taut=16), 96 (taut=10), 160 (taut=6) ...
% T0=taut;  %t_tau=Tmax/M; % tau=1/4 to 12
% nt=1;
% 
% %%---DHTSP setting ended ---

for run=1:runmax
    Archive=[];
    
    [start_weight,neighbour_indexMatrix]=init_weights(N,niche,objDim);

    g=1;
    Qop=init_Pop(varDim,N,bounds);
    % evaluate initial population
    Pop= func(Qop, 0); 

    idealpoint=ones(1,objDim)*inf;  % ideal point initialisation
    idealpoint=min(idealpoint,min(Pop(:,varDim+1:end), [],1)); 
    
    changetime=0;

    tStart=tic;
    while g<=Tg
        t=prob.update_time(g, taut, nt, T0); % update time instant
        checkvalue=checkChange(func, Pop,objDim, t, 0.1); % 10% population for change detection
        disp(['Evolve at generation ',num2str(g), ', time t= ', num2str(t), '; ----->']);
        
        if checkvalue==1  % in the event of change
            changetime=changetime+1;
            
            [IGD,SP,HV]=Performance(Pop(:,varDim+1:end),prob,t-1);%%% calculate performance indicators

            IGDt(changetime)=IGD;
            SPt(changetime)=SP;
            HVt(changetime)=HV;
    
            Archive=[Archive;Pop];%%%%
            
            Pop=svr_response(Archive,bounds,objDim,N,changetime,func,t);%SVR-style change response      
                
            idealpoint=min(Pop(:,varDim+1:end), [], 1);
            
        end
           
        for i=1:length(Pop)
            % one-by-one individual evolution
            indPop=genetic_operate(neighbour_indexMatrix,Pop,i,F,CR,bounds,varDim);
            indPop=checkbound(indPop,bounds);
            indPop=func(indPop,t);
            
            % update ideal point
            idealpoint=min(idealpoint,indPop(varDim+1:end));
            
            %%%% update neighbourhood
            neighbourindex=neighbour_indexMatrix(i,:);
            %  weights that neighbouring members have
            neighbour_weight=start_weight(neighbourindex,:);  
            %  neighbouring members
            neighbour_member=Pop(neighbourindex,:);
            
            % chebysheff decomposiiton
            new_tevalue=subjective_te(neighbour_weight,indPop,idealpoint,varDim);
            
            old_tevalue=subjective_te(neighbour_weight,neighbour_member,idealpoint,varDim);

            for k=1:niche % update at most niche members by new individual
                if new_tevalue(k)<old_tevalue(k)
                    Pop(neighbourindex(k),:)=indPop;
                end
            end
            
        end
        
        g=g+1;
    end
    %%%%%
    tEnd=toc(tStart);
    
    Time(run)=tEnd;
    IGDs(run,:)=IGDt;
    SPs(run,:)=SPt;
    HVs(run,:)=HVt;
    Archives{run}=Archive;

end
end


