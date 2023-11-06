function [Archives,Time, IGDs, HVs, SPs]=MOEADSVR(probID,N,runmax)
% MOEA/D-SVR for Dynamic Multi-objective Optimization
%
% Modified by Dr. Shouyong Jiang on 4/1/2023

%%%%%%%%Algorithmic parameters%%%%%%%%%%%%%%%%%%%%
niche=floor(N/5);%�����ģ
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
    %%%%%%%%%%%�����ʼȨ�ؾ���ÿ�������Ӧ���������ָ��
    [start_weight,neighbour_indexMatrix]=init_weights(N,niche,objDim);
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%��Ⱥ��������
    g=1;
    Qop=init_Pop(varDim,N,bounds);
    %%%%%%%�����ʼ��Ⱥ��Ӧ�ĺ���ֵ �����������objDimά
    Pop= func(Qop, 0); 
    %���������;
    
    idealpoint=ones(1,objDim)*inf;  %��ʼ��������Ϊ����
    idealpoint=min(idealpoint,min(Pop(:,varDim+1:end), [],1)); 
    
    changetime=0;

    tStart=tic;
    while g<=Tg
        t=prob.update_time(g, taut, nt, T0); % update time instant
        checkvalue=checkChange(func, Pop,objDim, t, 0.1); % 10% population for change detection
        disp(['Evolve at generation ',num2str(g), ', time t= ', num2str(t), '; ----->']);
        
        if checkvalue==1  %%˵���仯��
            changetime=changetime+1;
            
            [IGD,SP,HV]=Performance(Pop(:,varDim+1:end),prob,t-1);%%%����IGD��SP��

            IGDt(changetime)=IGD;
            SPt(changetime)=SP;
            HVt(changetime)=HV;

            Archive=[Archive;Pop];%%%%
            Pop=svr_response(Archive,bounds,objDim,N,changetime,func,t);%����Ӧ�Ա仯��������Ⱥ����        
                
            idealpoint=min(Pop(:,varDim+1:end), [], 1);
        end
           
        for i=1:length(Pop)
            %%%%ÿһ������ִ�н�����������һ���µĸ���
            indPop=genetic_operate(neighbour_indexMatrix,Pop,i,F,CR,bounds,varDim);
            indPop=checkbound(indPop,bounds);
            %��������
            indPop=func(indPop,t);
            
            %                 FES=FES+1;
            %%ÿ�β����¸�����������
            idealpoint=min(idealpoint,indPop(varDim+1:end));
            
            %%%%��������
            %��������ָ��
            neighbourindex=neighbour_indexMatrix(i,:);
            %             %%%�����Ӧ��Ȩ��
            neighbour_weight=start_weight(neighbourindex,:);  %�������ָ���Ӧ��Ȩ�ؾ���
            %             %%%�������
            neighbour_member=Pop(neighbourindex,:);
            
            %%%�¸����Ӧ���б�ѩ��ֵ
            
            new_tevalue=subjective_te(neighbour_weight,indPop,idealpoint,varDim);
%                         new_tevalue=subjective_ws(neighbour_weight,indPop,D,Dim);
            
            old_tevalue=subjective_te(neighbour_weight,neighbour_member,idealpoint,varDim);
%                         old_tevalue=subjective_ws(neighbour_weight,neighbour_member,D,Dim);
            for k=1:niche %���������
                if new_tevalue(k)<old_tevalue(k)
                    Pop(neighbourindex(k),:)=indPop;  %%%Ȼ��Ȩ�ز��䣬ֻ���滻�˸���
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


