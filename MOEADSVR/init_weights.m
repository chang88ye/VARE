function [subp,neighbour_indexMatrix]=init_weights(popsize,niche,ObjNum)
%��� Ȩ�ؾ��� �� ÿ�������Ӧ���������ָ��
subp=[];%%% Ȩ�ؾ���
neighbour_indexMatrix=[]; %%ÿ�������Ӧ���������ָ��
for i=1:popsize
    if ObjNum==2
        Weight=linspace(i,popsize-i,ObjNum);
        weight=Weight./popsize;
%         weight=zeros(1,2);
%         weight(1)=i/popsize;
%         weight(2)=(popsize-i)/popsize;
        subp=[subp;weight];
    elseif ObjNum==3
        Weight=linspace(i,popsize-i,ObjNum);
        weight=Weight./popsize;
        subp=[subp;weight];
    end
end

leng=length(subp);%�����нϴ��ָ��
distanceMatrix=zeros(leng,leng);

for i=1:leng
    for j=i+1:leng
        A=subp(i,:);
        B=subp(j,:);
        distanceMatrix(i,j)=(A-B)*(A-B)'; %�����������֮��ľ���
        distanceMatrix(j,i)= distanceMatrix(i,j);%���Գ�
    end
    [s,index]=sort(distanceMatrix(i,:));
    neighbour_indexMatrix(i,:)=index(1:niche);
end
        
         
        
        
        
        
        
        