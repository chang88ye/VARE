function Pop=pps_predit(TotalArchive,Popcenter,lu,ObjNum,t,popsize,Rch,func)
p=2;
M=23;
[~,SD]=size(TotalArchive);
D=SD-ObjNum;

if Rch<=p  %%%%˵���ǵ�һ�λ��ߵڶ���,һ����������һ���Ǽ̳��ϴ���Ⱥ����
    pop_current=TotalArchive((Rch-1)*popsize+1:Rch*popsize,:);%%��ǰ����Ⱥ
    Qop=pop_current(:,1:D);
    nr=round(popsize/2);
    Cop=init_Pop(D,nr,lu); %%%%��һ���������
    Sdex=randperm(popsize);
    Qop(Sdex(1:nr),1:D)=Cop;
else
    if Rch<=M
        seledex=1:M;
        %%��������ѡȡ
        selepop=Popcenter(seledex,:);

        for dd=1:D
            model=ar(selepop(:,dd),p,'fb');
            QZ=iddata(selepop(:,dd));%%%ֻ��23��
            Squar=model.Report.Fit.MSE;
            pp=predict(QZ(Rch-p:Rch),model,1);%%%%��3����Ԥ���4��QZֻ��23��
            newcenter(dd)=pp.OutputData(end)+normrnd(0,Squar);%%%1��D�ĸ��壬��������Ⱥ������
        end
    else
        seledex=(Rch-M):Rch;%%%%ǰ����Popcenter��ȺҪ���������������ٴ�1,��ѡһ��0
        selepop=Popcenter(seledex,:);
        for dd=1:D
            model=ar(selepop(:,dd),p,'fb');
            QZ=iddata(selepop(:,dd));%%%ֻ��23��
            Squar=model.Report.Fit.MSE;
            pp=predict(QZ(end-p:end),model,1);%%%%��3����Ԥ���4��QZֻ��23��
            newcenter(dd)=pp.OutputData(end)+normrnd(0,Squar);%%%1��D�ĸ���
        end
        %
    end

    %%%%%%%%%%%%%����Ԥ��,���α仯����Ⱥ����
    pop_before=TotalArchive((Rch-2)*popsize+1:(Rch-1)*popsize,:); %ǰһ����Ⱥ
    pop_current=TotalArchive((Rch-1)*popsize+1:Rch*popsize,:);%%��ǰ����Ⱥ

    manifold_current=pop_current-repmat(mean(pop_current,1),popsize,1);
    manifold_before=pop_before-repmat(mean(pop_before,1),popsize,1);
    cault_dis=pdist2(manifold_current(:,1:D),manifold_before(:,1:D));%^%%�����һ������
    [sdist,~]=sort(cault_dis,2);
    sumvalue=sum(sdist(:,1));%%%һ����ֵ
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

