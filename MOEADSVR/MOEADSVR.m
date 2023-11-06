function [Archives,Time, IGDs, HVs, SPs]=MOEADSVR(probID,N,runmax)
% MOEA/D-SVR for Dynamic Multi-objective Optimization
%
% Modified by Dr. Shouyong Jiang on 4/1/2023

%%%%%%%%Algorithmic parameters%%%%%%%%%%%%%%%%%%%%
niche=floor(N/5);%邻域规模
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

for run=1:runmax
    Archive=[];
    %%%%%%%%%%%输出初始权重矩阵及每个个体对应的邻域个体指标
    [start_weight,neighbour_indexMatrix]=init_weights(N,niche,objDim);
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%种群正常进化
    g=1;
    Qop=init_Pop(varDim,N,bounds);
    %%%%%%%计算初始种群对应的函数值 并保存在最后objDim维
    Pop= func(Qop, 0); 
    %更新理想点;
    
    idealpoint=ones(1,objDim)*inf;  %初始化理想点均为无穷
    idealpoint=min(idealpoint,min(Pop(:,varDim+1:end), [],1)); 
    
    changetime=0;

    tStart=tic;
    while g<=Tg
        t=prob.update_time(g, taut, nt, T0); % update time instant
        checkvalue=checkChange(func, Pop,objDim, t, 0.1); % 10% population for change detection
        disp(['Evolve at generation ',num2str(g), ', time t= ', num2str(t), '; ----->']);
        
        if checkvalue==1  %%说明变化了
            changetime=changetime+1;
            
            [IGD,SP,HV]=Performance(Pop(:,varDim+1:end),prob,t-1);%%%计算IGD与SP，

            IGDt(changetime)=IGD;
            SPt(changetime)=SP;
            HVt(changetime)=HV;

            Archive=[Archive;Pop];%%%%
            Pop=svr_response(Archive,bounds,objDim,N,changetime,func,t);%用于应对变化产生新种群个体        
                
            idealpoint=min(Pop(:,varDim+1:end), [], 1);
        end
           
        for i=1:length(Pop)
            %%%%每一个个体执行进化操作生成一个新的个体
            indPop=genetic_operate(neighbour_indexMatrix,Pop,i,F,CR,bounds,varDim);
            indPop=checkbound(indPop,bounds);
            %评估个体
            indPop=func(indPop,t);
            
            %                 FES=FES+1;
            %%每次产生新个体更新理想点
            idealpoint=min(idealpoint,indPop(varDim+1:end));
            
            %%%%更新邻域
            %个体邻域指标
            neighbourindex=neighbour_indexMatrix(i,:);
            %             %%%个体对应的权重
            neighbour_weight=start_weight(neighbourindex,:);  %邻域个体指标对应的权重矩阵
            %             %%%邻域个体
            neighbour_member=Pop(neighbourindex,:);
            
            %%%新个体对应的切贝雪夫值
            
            new_tevalue=subjective_te(neighbour_weight,indPop,idealpoint,varDim);
%                         new_tevalue=subjective_ws(neighbour_weight,indPop,D,Dim);
            
            old_tevalue=subjective_te(neighbour_weight,neighbour_member,idealpoint,varDim);
%                         old_tevalue=subjective_ws(neighbour_weight,neighbour_member,D,Dim);
            for k=1:niche %邻域个体数
                if new_tevalue(k)<old_tevalue(k)
                    Pop(neighbourindex(k),:)=indPop;  %%%然而权重不变，只是替换了个体
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


