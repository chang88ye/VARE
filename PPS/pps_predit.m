function Pop=pps_predit(TotalArchive,Popcenter,lu,ObjNum,t,popsize,Rch,func)
p=2;
M=23;
[~,SD]=size(TotalArchive);
D=SD-ObjNum;

if Rch<=p  %%%%说明是第一次或者第二次,一半是新生成一半是继承上代种群个体
    pop_current=TotalArchive((Rch-1)*popsize+1:Rch*popsize,:);%%当前代种群
    Qop=pop_current(:,1:D);
    nr=round(popsize/2);
    Cop=init_Pop(D,nr,lu); %%%%另一半随机产生
    Sdex=randperm(popsize);
    Qop(Sdex(1:nr),1:D)=Cop;
else
    if Rch<=M
        seledex=1:M;
        %%进行中心选取
        selepop=Popcenter(seledex,:);

        for dd=1:D
            model=ar(selepop(:,dd),p,'fb');
            QZ=iddata(selepop(:,dd));%%%只有23个
            Squar=model.Report.Fit.MSE;
            pp=predict(QZ(Rch-p:Rch),model,1);%%%%最3个来预测第4个QZ只有23个
            newcenter(dd)=pp.OutputData(end)+normrnd(0,Squar);%%%1乘D的个体，新生成种群的中心
        end
    else
        seledex=(Rch-M):Rch;%%%%前提是Popcenter种群要比最大迭代次数到少大1,多选一个0
        selepop=Popcenter(seledex,:);
        for dd=1:D
            model=ar(selepop(:,dd),p,'fb');
            QZ=iddata(selepop(:,dd));%%%只有23个
            Squar=model.Report.Fit.MSE;
            pp=predict(QZ(end-p:end),model,1);%%%%最3个来预测第4个QZ只有23个
            newcenter(dd)=pp.OutputData(end)+normrnd(0,Squar);%%%1乘D的个体
        end
        %
    end

    %%%%%%%%%%%%%流行预测,两次变化的种群即可
    pop_before=TotalArchive((Rch-2)*popsize+1:(Rch-1)*popsize,:); %前一代种群
    pop_current=TotalArchive((Rch-1)*popsize+1:Rch*popsize,:);%%当前代种群

    manifold_current=pop_current-repmat(mean(pop_current,1),popsize,1);
    manifold_before=pop_before-repmat(mean(pop_before,1),popsize,1);
    cault_dis=pdist2(manifold_current(:,1:D),manifold_before(:,1:D));%^%%计算出一个矩阵
    [sdist,~]=sort(cault_dis,2);
    sumvalue=sum(sdist(:,1));%%%一个数值
    Dvalue=sumvalue/(size(manifold_current,1));
    mdeta=Dvalue/D;%%Dvalue^2/D

    for kkk=1:D
        Qop(:,kkk)=manifold_current(:,kkk)+newcenter(kkk)+normrnd(0,mdeta);
    end
end

% bound-check and evaluate the new population
Qop=checkbound(Qop,lu);
Pop= func(Qop, t);
end

