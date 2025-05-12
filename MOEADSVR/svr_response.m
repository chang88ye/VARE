function Pop=svr_response(TotalArchive,lu,ObjNum,popsize,changetime,func,t)
q=4; % number of lags
pop_current=TotalArchive(end-popsize+1:end,:);% current population

[~,SD]=size(TotalArchive);
D=SD-ObjNum;

if changetime>q+1
    Qop=zeros(popsize, D);
    for i=1:popsize
        var=zeros(1,D);
        try
            for j=1:D
                % prepare training data for each variable of individual along the i-th weight vector
                train=[];
                for c =0:changetime-q-1
                    idx=popsize*((0:q)+c)+i;
                    train(end+1,:)=TotalArchive(idx,j);
                end
                x_train=train(:,1:end-1);
                y_train=train(:,end);

                % create a SVR model with Gausian RBF kernel function
                svr=fitrsvm(x_train, y_train, 'KernelFunction','Gaussian');%, 'OptimizeHyperparameters','auto'); %  (epsilon=0.05, C=1e3 in python implementation)
                var(j)=predict(svr,train(end, 2:end));
            end
            Qop(i,:)=var;
        catch
            Qop(i,:)=pop_current(i,1:D);
        end
        
    end
else
    Qop=pop_current(:,1:D);
end

% bound check and fitness evaluation
Qop=checkbound(Qop,lu);
Pop=func(Qop,t);
end
