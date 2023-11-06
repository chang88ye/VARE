function [Archives,Time, IGDs, HVs, SPs]=PPS(probID,N,runMax)

prob=problem(probID);
[varDim, objDim, bounds, func]=deal(prob.varDim, prob.objDim, prob.bounds, prob.func); % number of variables/objectives, variable bounds, and function handle of the problem

% get change severity/frequency and the past number of generations at which first change occurs

% method 1: get default values from problem definition
%[T0, taut, nt] =deal(instance.dynPar(1),instance.dynPar(2),instance.dynPar(3));

% method 2: get global values initialised in the main function
[nt, taut, T0]=dynamic_param(); 

Tg=T0+10*nt*taut; %set the termination genenerations: {T0 +10*nt*taut}

for run=1:runMax
    
    Archive=[];
    g=1;
    
    %%%%%%%%%%%%%%%%%% population initialisation
    Qop=init_population(prob,N);
    
    %%%%%%% evaluate initial population
    Pop=func(Qop,0);
    
    changetime=0;
    
    centerpop=zeros(Tg,varDim);

    %%%%%%%%%%%%%%%%%%%%%Optimization
    tStart=tic;
    while g<=Tg

        t=prob.update_time(g, taut, nt, T0); % update time instant
        disp(['Evolve at generation ',num2str(g), ', time t= ', num2str(t), '; ----->']);

        checkvalue=checkChange(func, Pop,objDim, t, 0.1); % 10% population for change detection
        
        %%%%%%%%%%%%%%%%%%%%%%Response_change%%%%%%%%%%%%%%%%
        if checkvalue==1  
            changetime=changetime+1;
            t_prev=prob.update_time(g-1, taut, nt, T0);

            [IGD,SP,HV]= Performance(Pop(:, end-objDim+1:end),prob,t_prev);
            
            IGDt(changetime)=IGD;
            SPt(changetime)=SP;
            HVt(changetime)=HV;

            Archive=[Archive;Pop];%%%%
            
            %%%%% saving centroids in an extenal archive
            centerpop(changetime,:)=mean(Pop(:,1:varDim),1);
            
            Qop=pps_predit(Archive,centerpop,bounds,objDim,t,N,changetime,func);
            

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%  Evolutionary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        QWE = RMMEAD_operator(Pop,varDim,objDim);
        Qop=checkbound(QWE,bounds);
        Offspring= func(Qop, t); %% evaluate the new population

        %%%%%%%%%%%%%%%%%%%%%%   selection the offspring from mix pop
        Combin_Population=[Pop;Offspring];
        Zop = EnvironmentalSelection(Combin_Population,N,varDim,objDim);
        Pop=Zop;
        
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


