function [W, K, n1, n2]=getW(Xs, Xt, mu, lambda, dim, kind, p1, p2, p3)
%ԭTCA����
%Xs��Դ������
%Xt��Ŀ�������ݣ���Xs������ͬ
%mu��ƽ�����ӣ�Խ��Խ����ӳ�������ƶȣ�ԽСԽ����W�ĸ��Ӷ�
%lambda����������ʹ�øò�����ֻ��Ϊ����ʽ��ͳһ
%dim����dimΪ���ڵ���1������ʱ��dimΪ��ά��Ŀ��ά����
%     ��dimΪ����0С��1��С��ʱ����ȡ����������Ӧ������ֵ�ĺ�>=ȫ������ֵ�Ӻ�*dim
%kind���˺���ѡ��:'Gaussian'��'Laplacian'��'Polynomial',����һ�ɷ���-1
%p1,p2,p3���˺�����Ҫ�����Ĳ���
%W���任����n1+n2->dim
%K�����任����
%n1,n2��Դ���ݣ�Ŀ�����ݵ���Ŀ

    n1 = size(Xs, 2);
	n2 = size(Xt, 2);
    
%%%%%%%%%%% ����K %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X(:,1:n1)=Xs;
    X(:,n1+1:n1+n2)=Xt;
    for i=1:n1+n2 
        for j=1:n1+n2 
            K(i,j)=getKernel(X(:,i), X(:,j), kind, p1, p2, p3);
        end
    end
    
%%%%%%%%%%% ����L %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    L(1:n1, 1:n1)=ones(n1, n1)/(n1*n1);
    L(n1+1:n1+n2, n1+1:n1+n2)=ones(n2, n2)/(n2*n2);
    L(1:n1, n1+1:n1+n2)=ones(n1, n2)/(-n1*n2);
    L(n1+1:n1+n2, 1:n1)=ones(n2, n1)/(-n1*n2);
    
%%%%%%%%%%% ����H %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    H=eye(n1+n2)-ones(n1+n2, n1+n2)/(n1+n2);
    
%%%%%%%%%%% ����W %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Temp=(eye(n1+n2)+mu*K*L*K)^(-1)*K*H*K;
    [V,D]=eig(Temp);
    V = real(V);
    D = real(D);
    D=diag(D);
    [D,I]=sort(D,'descend');
    
    if dim>0 && dim<1
        count=0;
        cur=0;
        s=sum(D);
        while cur/s<dim
            count=count+1;
            cur=cur+D(count);
        end
    else
        count=dim;
    end
    
    I=I(1:count,1);
    W=V(:,I');
end

    