% cec09.m
% 
% The Matlab version of the test instances for CEC 2009 Multiobjective
%   Optimization Competition.
% 
% Usage: fobj = cec09(problem_name), the handle of the function will be
%   with fobj
% 
% Please refer to the report for correct one if the source codes are not
%   consist with the report.
% History:
%   v1 Sept.08 2008
%   v2 Nov.18  2008
%   v3 Nov.26  2008

function fobj = cec09(name)
 %%%%%约束条件的写法均为大于0
    switch name
        case 1
            fobj = @DHTPS;%%%%Dynamic Hydro-Thermal Power Scheduling Problem (DHPSP)
        case 2
            fobj = @DSRD; %%Dynamic Speed reducer design problem
        case 3
            fobj = @DWBD;%Welded beam design problem 
        case 4
            fobj = @DBD;%%%%%Disk brake design problem 
        case 5
            fobj = @CBD;%%% Cantilever beam design problem
        case 6
            fobj = @UF6;
        case 7
            fobj = @UF7;
        case 8
            fobj = @UF8; 
        case 9
            fobj = @UF9; 
        case 10
            fobj = @UF10;
        case 11
            fobj = @CF1;
        case 12
            fobj = @CF2; 
        case 13
            fobj = @CF3;  
        case 14
            fobj = @CF4;
        case 15
            fobj = @CF5; 
        case 16
            fobj = @CF6;
        case 17
            fobj = @CF7;
        case 18
            fobj = @CF8; 
        case 19
            fobj = @CF9; 
        case 20
            fobj = @CF10;           
        otherwise
            fobj = @UF1;
    end
end

%% Dynamic Hydro-Thermal Power Scheduling Problem (DHPSP)
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function [y,c] = DHTPS(x,t)
    %%%%有变量次数t有关的控制参数
    if t==0
        Pdt=900;
    elseif t==1
        Pdt=1000;
    elseif t==2
        Pdt=1100;
    else
        Pdt=1300;
    end
        
    [dim, num]   = size(x);
    y            = zeros(dim,num);
    
    k3=60+ 1.8*x(3,:)+0.003*x(3,:).^2 +abs(140*sin(0.04 *20*x(3,:)));
    k4=100+2.1*x(4,:)+0.0012*x(4,:).^2+abs(160*sin(0.038*30*x(4,:)));
    k5=120+2.0*x(5,:)+0.0010*x(5,:).^2+abs(180*sin(0.037*40*x(5,:)));
    k6=40 +1.8*x(6,:)+0.0015*x(6,:).^2+abs(200*sin(0.035*50*x(6,:)));    
    y(1,:)       = 12*(k3+k4+k5+k6);
    
    kk3=50-0.555*x(3,:)+0.0150*x(3,:).^2+0.5773*exp(0.02446*x(3,:));
    kk4=60-1.355*x(4,:)+0.0105*x(4,:).^2+0.4968*exp(0.02270*x(4,:));
    kk5=45-0.600*x(5,:)+0.0080*x(5,:).^2+0.4860*exp(0.01948*x(5,:));
    kk6=30-0.555*x(6,:)+0.0120*x(6,:).^2+0.5035*exp(0.02075*x(6,:));
    
    y(2,:)       = 12*(kk3+kk4+kk5+kk6);   
    
    B=10^(-6)*[49 14 15 15 20 17;14 45 16 20 18 15;15 16 39 10 12 12;15 20 10 40 14 10;20 18 12 14 35 11; 17 15 12 10 11 36];
    Ph=0;
   for z=1:6
        PPh=0;
        for zz=1:6
        PPh=PPh+x(z,:)*B(z,zz).*x(zz,:);
        end
        Ph=Ph+PPh;
   end

    c(1,:)       = x(3,:)+x(4,:)+x(5,:)+x(6,:)+x(1,:)+x(2,:)-Pdt-Ph;
    c(2,:)       = 12*(260+8.5*x(1,:)+0.00986*x(1,:).^2)-125/48;   
    c(3,:)       = 12*(250+9.8*x(2,:)+0.01140*x(2,:).^2)-286/48;   
    
end
%% Speed reducer design problem
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function [y,c] = DSRD(x,t)
    %%%%%均与t与关的控制参数 
    a1=0.7854+t/10;
    a2=14.933+t/10;
    a3=43.0934+t/10;
    a4=1.508+t/10;
    a5=7.477+t/10;
    a6=0.7854+t/10;
    a7=0.1+0.25-1/(t+4);
    
    [dim, num]   = size(x);
    Y            = zeros(dim,num);

    p1=10*(x(3,:).^2)/3+a2*x(3,:)-a3;
    p2=x(6,:).^2+x(7,:).^2;
    p3=x(6,:).^3+x(7,:).^3;
    p4=x(4,:).*(x(6,:).^2)+x(5,:).*(x(7,:).^2);
    y(1,:)       = (a1*x(1,:).*(x(2,:).^2).*p1-a4*x(1,:).*p2+a5*p3+a6*p4);
    
    p5=745*x(4,:)./(x(2,:).*x(3,:));
    y(2,:)       = sqrt(p5.^2+1.69*(10^7))./(a7*x(6,:).^3);
    
    
    %%%约束条件大于0
    c(1,:)       = 1/27-1./(x(1,:).*(x(2,:).^2).*x(3,:));
    c(2,:)       = 1/397.5-1./(x(1,:).*(x(2,:).^2).*(x(3,:).^2));
    c(3,:)       = 1/1.93-(x(4,:).^3)./(x(2,:).*x(3,:).*(x(6,:).^4));
    c(4,:)       = 1/1.93-(x(5,:).^3)./(x(2,:).*x(3,:).*(x(7,:).^4));
    c(5,:)       = 40-x(2,:).*x(3,:);
    c(6,:)       = 12-x(1,:)./x(2,:);
    c(7,:)       = x(1,:)./x(2,:)-5;
    c(8,:)       = x(4,:)-1.9-1.5*x(6,:);
    c(9,:)       = x(5,:)-1.9-1.1*x(1,:);
    c(10,:)       = 4300-y(1,:);
    p6=745*x(5,:)./(x(2,:).*x(3,:));
    c(11,:)       = 1100-sqrt(p6.^2+1.575*(10^8))/(0.1*x(7,:).^3);
    clear Y;
end
%%%Welded beam design problem  
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function [y,c] = DWBD(x,t)
    %%%%有变量次数t有关的控制参数
    if t==0
        Pt=10000;
    elseif t==1
        Pt=8000;
    elseif t==2
        Pt=6000;
    else
        Pt=3000;
    end
        
    [dim, num]   = size(x);
    Y            = zeros(dim,num);
    
    y(1,:)       = (1+0.10471)*(x(1,:).^2).*x(2,:)+0.04811*x(3,:).*x(4,:).*(14+x(2,:));
    y(2,:)       = (4*Pt*14^3)/(3*(10^7)*x(3,:).^3.*x(4,:));   
    
    
    tau1=Pt./(sqrt(2)*x(1,:).*x(2,:));
    
    M=Pt*(14+0.5*x(2,:));
    R=sqrt((x(2,:)^2/4)+((x(1,:)+x(3,:))/2)^2);
    J=2*sqrt(2)*x(1,:).*x(2,:).*((x(2,:).^2/12)+((x(1,:)+x(3,:))/2)^2);
    tau2=(M.*R)./J;
    
    tau=sqrt(tau1.^2+2*tau1.*tau2.*(0.5*x(2,:)./R)+tau2.^2);
      
    deltax=6*Pt*14./(x(3,:).^2.*x(4,:));
    
    pcx1=(4.013*1.2*10^7*sqrt(x(3,:).^2.*x(4,:).^6/36))/(14^2);
    pcx2=(x(3,:)./28).*(sqrt(3*10^7/(4*1.2*10^7)));
    pcx=pcx1.*(1-pcx2);
    
    c(1,:)       = 13600-tau;
    c(2,:)       = 30000-deltax;   
    c(3,:)       = x(4,:)-x(1,:);   
    c(4,:)       = pcx-Pt;
     
    clear Y;
end

%% UF4
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function [y,c] = DBD(x)
    
    [dim, num]   = size(x);
    Y            = zeros(dim,num);
    p1=x(2,:).^2-x(1,:).^2;
    y(1,:)       = 4.9*(10^(-5))*p1.*(x(4,:)-1);
    p2=x(2,:).^3-x(1,:).^3;    
    y(2,:)       = 9.82*(10^6)*p1./(x(3,:).*x(4,:).*p2);
    
    c(1,:)       = x(2,:)-x(1,:)-20;
    c(2,:)       = 30-2.5*(x(4,:)+1);   
    c(3,:)       = 0.4-x(3,:)./(3.14*p1.^2);   
    c(4,:)       = 1-(2.22*(10^(-3))*x(3,:).*p2)./(p1.^2);
    c(5,:)       = (2.66*(10^(-2))*x(3,:).*x(4,:).*p2)./(p1)-900;    
    clear Y;
end


%% Cantilever beam design problem这个还是有问题
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function [y,c] = CBD(x)
    
    [dim, num]   = size(x);
    Y            = zeros(dim,num);
    
    y(1,:)       = 0.25*7800*pi*x(2,:).*(x(1,:).^2);
    p1=64*(x(2,:).^3);
    p2=3*207000000*pi*(x(1,:).^4);
    y(2,:)       = p1./p2;
    
    c(1,:)       = 30000-((32*x(2,:))./(pi*(x(1,:).^3)));
    c(2,:)       = 0.005-p1./p2;   
       
    clear Y;
end
%% Four-bar truss design problem
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function y = FBT(x)
    [dim, num]  = size(x);
    tmp         = zeros(dim,num);
    y(1,:)      = 200*(2*x(1,:)+sqrt(2*x(2,:))+sqrt(x(2,:))+x(4,:));
    p1=2./x(1,:);
    p2=2*sqrt(2)./x(2,:);
    p3=2*sqrt(2)./x(3,:);
    p4=2./x(4,:);
    y(2,:)      = 0.01*(p1+p2-p3+p4);
    clear tmp;
end

%% UF6
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function y = UF6(x)
    N            = 2.0;
    E            = 0.1;
    [dim, num]   = size(x);
    Y            = zeros(dim,num);
    Y(2:dim,:)  = x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));
    tmp1         = zeros(dim,num);
    tmp1(2:dim,:)= Y(2:dim,:).^2;
    tmp2         = zeros(dim,num);
    tmp2(2:dim,:)= cos(20.0*pi*Y(2:dim,:)./sqrt(repmat((2:dim)',[1,num])));
    tmp11        = 4.0*sum(tmp1(3:2:dim,:)) - 2.0*prod(tmp2(3:2:dim,:)) + 2.0;  % odd index
    tmp21        = 4.0*sum(tmp1(2:2:dim,:)) - 2.0*prod(tmp2(2:2:dim,:)) + 2.0;  % even index
    tmp          = max(0,(1.0/N+2.0*E)*sin(2.0*N*pi*x(1,:)));
    y(1,:)       = x(1,:)       + tmp + 2.0*tmp11/size(3:2:dim,2);
    y(2,:)       = 1.0 - x(1,:) + tmp + 2.0*tmp21/size(2:2:dim,2);
    clear Y tmp1 tmp2;
end

%% UF7
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function y = UF7(x)
    [dim, num]  = size(x);
    Y           = zeros(dim,num);
    Y(2:dim,:)  = (x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]))).^2;
    tmp1        = sum(Y(3:2:dim,:));  % odd index
    tmp2        = sum(Y(2:2:dim,:));  % even index
    tmp         = (x(1,:)).^0.2;
    y(1,:)      = tmp       + 2.0*tmp1/size(3:2:dim,2);
    y(2,:)      = 1.0 - tmp + 2.0*tmp2/size(2:2:dim,2);
    clear Y;
end

%% UF8
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function y = UF8(x)
    [dim, num]  = size(x);
    Y           = zeros(dim,num);
    Y(3:dim,:)  = (x(3:dim,:) - 2.0*repmat(x(2,:),[dim-2,1]).*sin(2.0*pi*repmat(x(1,:),[dim-2,1]) + pi/dim*repmat((3:dim)',[1,num]))).^2;
    tmp1        = sum(Y(4:3:dim,:));  % j-1 = 3*k
    tmp2        = sum(Y(5:3:dim,:));  % j-2 = 3*k
    tmp3        = sum(Y(3:3:dim,:));  % j-0 = 3*k
    y(1,:)      = cos(0.5*pi*x(1,:)).*cos(0.5*pi*x(2,:)) + 2.0*tmp1/size(4:3:dim,2);
    y(2,:)      = cos(0.5*pi*x(1,:)).*sin(0.5*pi*x(2,:)) + 2.0*tmp2/size(5:3:dim,2);
    y(3,:)      = sin(0.5*pi*x(1,:))                     + 2.0*tmp3/size(3:3:dim,2);
    clear Y;
end

%% UF9
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function y = UF9(x)
    E           = 0.1;
    [dim, num]  = size(x);
    Y           = zeros(dim,num);
    Y(3:dim,:)  = (x(3:dim,:) - 2.0*repmat(x(2,:),[dim-2,1]).*sin(2.0*pi*repmat(x(1,:),[dim-2,1]) + pi/dim*repmat((3:dim)',[1,num]))).^2;
    tmp1        = sum(Y(4:3:dim,:));  % j-1 = 3*k
    tmp2        = sum(Y(5:3:dim,:));  % j-2 = 3*k
    tmp3        = sum(Y(3:3:dim,:));  % j-0 = 3*k
    tmp         = max(0,(1.0+E)*(1-4.0*(2.0*x(1,:)-1).^2));
    y(1,:)      = 0.5*(tmp+2*x(1,:)).*x(2,:)     + 2.0*tmp1/size(4:3:dim,2);
    y(2,:)      = 0.5*(tmp-2*x(1,:)+2.0).*x(2,:) + 2.0*tmp2/size(5:3:dim,2);
    y(3,:)      = 1-x(2,:)                       + 2.0*tmp3/size(3:3:dim,2);
    clear Y;
end

%% UF10
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function y = UF10(x)
    [dim, num]  = size(x);
    Y           = zeros(dim,num);
    Y(3:dim,:)  = x(3:dim,:) - 2.0*repmat(x(2,:),[dim-2,1]).*sin(2.0*pi*repmat(x(1,:),[dim-2,1]) + pi/dim*repmat((3:dim)',[1,num]));
    H           = zeros(dim,num);
    H(3:dim,:)  = 4.0*Y(3:dim,:).^2 - cos(8.0*pi*Y(3:dim,:)) + 1.0;
    tmp1        = sum(H(4:3:dim,:));  % j-1 = 3*k
    tmp2        = sum(H(5:3:dim,:));  % j-2 = 3*k
    tmp3        = sum(H(3:3:dim,:));  % j-0 = 3*k
    y(1,:)      = cos(0.5*pi*x(1,:)).*cos(0.5*pi*x(2,:)) + 2.0*tmp1/size(4:3:dim,2);
    y(2,:)      = cos(0.5*pi*x(1,:)).*sin(0.5*pi*x(2,:)) + 2.0*tmp2/size(5:3:dim,2);
    y(3,:)      = sin(0.5*pi*x(1,:))                     + 2.0*tmp3/size(3:3:dim,2);
    clear Y H;
end

%% CF1
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function [y,c] = CF1(x)
    a            = 1.0;
    N            = 10.0;
    [dim, num]   = size(x);
    Y            = zeros(dim,num);
    Y(2:dim,:)   = (x(2:dim,:) - repmat(x(1,:),[dim-1,1]).^(0.5+1.5*(repmat((2:dim)',[1,num])-2.0)/(dim-2.0))).^2;
    tmp1         = sum(Y(3:2:dim,:));% odd index
    tmp2         = sum(Y(2:2:dim,:));% even index 
    y(1,:)       = x(1,:)       + 2.0*tmp1/size(3:2:dim,2);
    y(2,:)       = 1.0 - x(1,:) + 2.0*tmp2/size(2:2:dim,2);
    c(1,:)       = y(1,:) + y(2,:) - a*abs(sin(N*pi*(y(1,:)-y(2,:)+1.0))) - 1.0;
    clear Y;
end

%% CF2
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function [y,c] = CF2(x)
    a           = 1.0;
    N           = 2.0;
    [dim, num]  = size(x);
    tmp         = zeros(dim,num);
    tmp(2:dim,:)= (x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]))).^2;
    tmp1        = sum(tmp(3:2:dim,:));  % odd index
    tmp(2:dim,:)= (x(2:dim,:) - cos(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]))).^2;
    tmp2        = sum(tmp(2:2:dim,:));  % even index
    y(1,:)      = x(1,:)             + 2.0*tmp1/size(3:2:dim,2);
    y(2,:)      = 1.0 - sqrt(x(1,:)) + 2.0*tmp2/size(2:2:dim,2);
    t           = y(2,:) + sqrt(y(1,:)) - a*sin(N*pi*(sqrt(y(1,:))-y(2,:)+1.0)) - 1.0;
    c(1,:)      = sign(t).*abs(t)./(1.0+exp(4.0*abs(t)));
    clear tmp;
end

%% CF3
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function [y,c] = CF3(x)
    a            = 1.0;
    N            = 2.0;
    [dim, num]   = size(x);
    Y            = zeros(dim,num);
    Y(2:dim,:)   = x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));
    tmp1         = zeros(dim,num);
    tmp1(2:dim,:)= Y(2:dim,:).^2;
    tmp2         = zeros(dim,num);
    tmp2(2:dim,:)= cos(20.0*pi*Y(2:dim,:)./sqrt(repmat((2:dim)',[1,num])));
    tmp11        = 4.0*sum(tmp1(3:2:dim,:)) - 2.0*prod(tmp2(3:2:dim,:)) + 2.0;  % odd index
    tmp21        = 4.0*sum(tmp1(2:2:dim,:)) - 2.0*prod(tmp2(2:2:dim,:)) + 2.0;  % even index
    y(1,:)       = x(1,:)          + 2.0*tmp11/size(3:2:dim,2);
    y(2,:)       = 1.0 - x(1,:).^2 + 2.0*tmp21/size(2:2:dim,2);
    c(1,:)       = y(2,:) + y(1,:).^2 - a*sin(N*pi*(y(1,:).^2-y(2,:)+1.0)) - 1.0;   
    clear Y tmp1 tmp2;
end

%% CF4
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function [y,c] = CF4(x)
    [dim, num]  = size(x);
    tmp         = zeros(dim,num);
    tmp(2:dim,:)= x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));
    tmp1        = sum(tmp(3:2:dim,:).^2);  % odd index
    tmp2        = sum(tmp(4:2:dim,:).^2);  % even index
    index1      = tmp(2,:) < (1.5-0.75*sqrt(2.0));
    index2      = tmp(2,:)>= (1.5-0.75*sqrt(2.0));
    tmp(2,index1) = abs(tmp(2,index1));
    tmp(2,index2) = 0.125 + (tmp(2,index2)-1.0).^2;
    y(1,:)      = x(1,:)                  + tmp1;
    y(2,:)      = 1.0 - x(1,:) + tmp(2,:) + tmp2;
    t           = x(2,:) - sin(6.0*pi*x(1,:)+2.0*pi/dim) - 0.5*x(1,:) + 0.25;
    c(1,:)      = sign(t).*abs(t)./(1.0+exp(4.0*abs(t)));
    clear tmp index1 index2;
end

%% CF5
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function [y,c] = CF5(x)
    [dim, num]  = size(x);
    tmp         = zeros(dim,num);
    tmp(2:dim,:)= x(2:dim,:) - 0.8*repmat(x(1,:),[dim-1,1]).*cos(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));
    tmp1        = sum(2.0*tmp(3:2:dim,:).^2-cos(4.0*pi*tmp(3:2:dim,:))+1.0);  % odd index
    tmp(2:dim,:)= x(2:dim,:) - 0.8*repmat(x(1,:),[dim-1,1]).*sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));    
    tmp2        = sum(2.0*tmp(4:2:dim,:).^2-cos(4.0*pi*tmp(4:2:dim,:))+1.0);  % even index
    index1      = tmp(2,:) < (1.5-0.75*sqrt(2.0));
    index2      = tmp(2,:)>= (1.5-0.75*sqrt(2.0));
    tmp(2,index1) = abs(tmp(2,index1));
    tmp(2,index2) = 0.125 + (tmp(2,index2)-1.0).^2;
    y(1,:)      = x(1,:)                  + tmp1;
    y(2,:)      = 1.0 - x(1,:) + tmp(2,:) + tmp2;
    c(1,:)      = x(2,:) - 0.8*x(1,:).*sin(6.0*pi*x(1,:)+2.0*pi/dim) - 0.5*x(1,:) + 0.25;
    clear tmp;
end

%% CF6
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function [y,c] = CF6(x)
    [dim, num]  = size(x);
    tmp         = zeros(dim,num);
    tmp(2:dim,:)= x(2:dim,:) - 0.8*repmat(x(1,:),[dim-1,1]).*cos(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));
    tmp1        = sum(tmp(3:2:dim,:).^2);  % odd index
    tmp(2:dim,:)= x(2:dim,:) - 0.8*repmat(x(1,:),[dim-1,1]).*sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));    
    tmp2        = sum(tmp(2:2:dim,:).^2);  % even index
    y(1,:)      = x(1,:)            + tmp1;
    y(2,:)      = (1.0 - x(1,:)).^2 + tmp2;
    tmp         = 0.5*(1-x(1,:))-(1-x(1,:)).^2;
    c(1,:)      = x(2,:) - 0.8*x(1,:).*sin(6.0*pi*x(1,:)+2*pi/dim) - sign(tmp).*sqrt(abs(tmp));
    tmp         = 0.25*sqrt(1-x(1,:))-0.5*(1-x(1,:));
    c(2,:)      = x(4,:) - 0.8*x(1,:).*sin(6.0*pi*x(1,:)+4*pi/dim) - sign(tmp).*sqrt(abs(tmp));    
    clear tmp;
end

%% CF7
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function [y,c] = CF7(x)
    [dim, num]  = size(x);
    tmp         = zeros(dim,num);
    tmp(2:dim,:)= x(2:dim,:) - cos(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));
    tmp1        = sum(2.0*tmp(3:2:dim,:).^2-cos(4.0*pi*tmp(3:2:dim,:))+1.0);  % odd index
    tmp(2:dim,:)= x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));
    tmp2        = sum(2.0*tmp(6:2:dim,:).^2-cos(4.0*pi*tmp(6:2:dim,:))+1.0);  % even index
    tmp(2,:)    = tmp(2,:).^2;
    tmp(4,:)    = tmp(4,:).^2;
    y(1,:)      = x(1,:)                                  + tmp1;
    y(2,:)      = (1.0 - x(1,:)).^2 + tmp(2,:) + tmp(4,:) + tmp2;
    tmp         = 0.5*(1-x(1,:))-(1-x(1,:)).^2;
    c(1,:)      = x(2,:) - sin(6.0*pi*x(1,:)+2*pi/dim) - sign(tmp).*sqrt(abs(tmp));
    tmp         = 0.25*sqrt(1-x(1,:))-0.5*(1-x(1,:));
    c(2,:)      = x(4,:) - sin(6.0*pi*x(1,:)+4*pi/dim) - sign(tmp).*sqrt(abs(tmp));    
    clear tmp;
end

%% CF8
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function [y,c] = CF8(x)
    N           = 2.0;
    a           = 4.0;
    [dim, num]  = size(x);
    Y           = zeros(dim,num);
    Y(3:dim,:)  = (x(3:dim,:) - 2.0*repmat(x(2,:),[dim-2,1]).*sin(2.0*pi*repmat(x(1,:),[dim-2,1]) + pi/dim*repmat((3:dim)',[1,num]))).^2;
    tmp1        = sum(Y(4:3:dim,:));  % j-1 = 3*k
    tmp2        = sum(Y(5:3:dim,:));  % j-2 = 3*k
    tmp3        = sum(Y(3:3:dim,:));  % j-0 = 3*k
    y(1,:)      = cos(0.5*pi*x(1,:)).*cos(0.5*pi*x(2,:)) + 2.0*tmp1/size(4:3:dim,2);
    y(2,:)      = cos(0.5*pi*x(1,:)).*sin(0.5*pi*x(2,:)) + 2.0*tmp2/size(5:3:dim,2);
    y(3,:)      = sin(0.5*pi*x(1,:))                     + 2.0*tmp3/size(3:3:dim,2);
    c(1,:)      = (y(1,:).^2+y(2,:).^2)./(1.0-y(3,:).^2) - a*abs(sin(N*pi*((y(1,:).^2-y(2,:).^2)./(1.0-y(3,:).^2)+1.0))) - 1.0;
    clear Y;
end

%% CF9
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function [y,c] = CF9(x)
    N           = 2.0;
    a           = 3.0;
    [dim, num]  = size(x);
    Y           = zeros(dim,num);
    Y(3:dim,:)  = (x(3:dim,:) - 2.0*repmat(x(2,:),[dim-2,1]).*sin(2.0*pi*repmat(x(1,:),[dim-2,1]) + pi/dim*repmat((3:dim)',[1,num]))).^2;
    tmp1        = sum(Y(4:3:dim,:));  % j-1 = 3*k
    tmp2        = sum(Y(5:3:dim,:));  % j-2 = 3*k
    tmp3        = sum(Y(3:3:dim,:));  % j-0 = 3*k
    y(1,:)      = cos(0.5*pi*x(1,:)).*cos(0.5*pi*x(2,:)) + 2.0*tmp1/size(4:3:dim,2);
    y(2,:)      = cos(0.5*pi*x(1,:)).*sin(0.5*pi*x(2,:)) + 2.0*tmp2/size(5:3:dim,2);
    y(3,:)      = sin(0.5*pi*x(1,:))                     + 2.0*tmp3/size(3:3:dim,2);
    c(1,:)      = (y(1,:).^2+y(2,:).^2)./(1.0-y(3,:).^2) - a*sin(N*pi*((y(1,:).^2-y(2,:).^2)./(1.0-y(3,:).^2)+1.0)) - 1.0;
    clear Y;
end

%% CF10
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function [y,c] = CF10(x)
    a           = 1.0;
    N           = 2.0;
    [dim, num]  = size(x);
    Y           = zeros(dim,num);
    Y(3:dim,:)  = x(3:dim,:) - 2.0*repmat(x(2,:),[dim-2,1]).*sin(2.0*pi*repmat(x(1,:),[dim-2,1]) + pi/dim*repmat((3:dim)',[1,num]));
    H           = zeros(dim,num);
    H(3:dim,:)  = 4.0*Y(3:dim,:).^2 - cos(8.0*pi*Y(3:dim,:)) + 1.0;
    tmp1        = sum(H(4:3:dim,:));  % j-1 = 3*k
    tmp2        = sum(H(5:3:dim,:));  % j-2 = 3*k
    tmp3        = sum(H(3:3:dim,:));  % j-0 = 3*k
    y(1,:)      = cos(0.5*pi*x(1,:)).*cos(0.5*pi*x(2,:)) + 2.0*tmp1/size(4:3:dim,2);
    y(2,:)      = cos(0.5*pi*x(1,:)).*sin(0.5*pi*x(2,:)) + 2.0*tmp2/size(5:3:dim,2);
    y(3,:)      = sin(0.5*pi*x(1,:))                     + 2.0*tmp3/size(3:3:dim,2);
    c(1,:)      = (y(1,:).^2+y(2,:).^2)./(1.0-y(3,:).^2) - a*sin(N*pi*((y(1,:).^2-y(2,:).^2)./(1.0-y(3,:).^2)+1.0)) - 1.0;
    clear Y H;
end