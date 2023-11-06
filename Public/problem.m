function prob=problem(prob_id)
%% Implementation of popular dynamic test suites
% 
% DF test suite:
%               14 test functions are for cec2018 competition on  |
%               dynamic multiobjective optimisation.              |
% SDP test suite:
%               15 objective/variable-scalable test problems
%
% FDA      :  ID 201 - 205
%             - The 5 test functions FDA1 - FDA5 pioneered research on EDMO
% 
% DMOP     :  ID 301 - 303
%             - three test functions DMOP1 - DMOP3 in dCOEA
%           
% F        :  ID 401 - 408
%             - Eight test functions F1 - F8 in PPS
    
% -----------------------------------------------------------------------|%
% The "time" term in the test suite is defined as:                       |
%          t=1/nt*floor(tau/taut)                                        |
% where - nt:    severity of change                                      |
%       - taut:  frequency of change                                     |
%       - tau:   current generation counter                              |
% -----------------------------------------------------------------------|%
%
% Any questions can be directed to          
%    Dr. Shouyong Jiang at math4neu@gmail.com.
% the first change occurs after T0 generations, that is, the
%generation at which a change occurs is (T0+1), (T0+taut+1), etc. 
% This document is free to disseminate for academic use. 
% test suites are numbered, i.e. CEC2018-DF (1), FDA (2), dCOEA-DMOP (3), PPS-F (4), etc
%
% 

    prob=struct('dynPar', [100,  10,  10]); % create a struct and provide values for TO, taut, and nt
    
    prob.update_time=@cal_time;
    prob.curT=0; % initialise time instant

    switch prob_id

        %% DF test suite: CEC2018
        case 101                 % 'DF1' (dMOP2)
            prob.name='DF1';
            prob.func=@DF1;
            prob.varDim=10;
            prob.objDim=2;
            prob.bounds=[zeros(1,prob.varDim);ones(1,prob.varDim)];
        case 102                 % DF2 (modified dMOP3)
            prob.name='DF2';
            prob.func=@DF2;
            prob.varDim=10;
            prob.objDim=2;
            prob.bounds=[zeros(1,prob.varDim);ones(1,prob.varDim)];
        case 103                 % 'DF3' (ZJZ)
            prob.name='DF3';
            prob.func=@DF3;
            prob.varDim=10;
            prob.objDim=2;
            prob.bounds=[0 -ones(1,prob.varDim-1);1 2*ones(1,prob.varDim-1)];
        case 104                 % DF4
            prob.name='DF4';
            prob.func=@DF4;
            prob.varDim=10;
            prob.objDim=2;
            prob.bounds=[-2*ones(1,prob.varDim);2*ones(1,prob.varDim)];
        case 105                 % DF5 (modified JY2)
            prob.name='DF5';
            prob.func=@DF5;
            prob.varDim=10;
            prob.objDim=2;
            prob.bounds=[0 -ones(1,prob.varDim-1);1 ones(1,prob.varDim-1)];
        case 106                 % DF6  (modified JY7)
            prob.name='DF6';
            prob.func=@DF6;
            prob.varDim=10;
            prob.objDim=2;
            prob.bounds=[zeros(1,prob.varDim);ones(1,prob.varDim)];
        case 107                 % DF7
            prob.name='DF7';
            prob.func=@DF7;
            prob.varDim=10;
            prob.objDim=2;
            prob.bounds=[1 zeros(1,prob.varDim-1);4 ones(1,prob.varDim-1)];
        case 108                 % DF8
            prob.name='DF8';
            prob.func=@DF8;
            prob.varDim=10;
            prob.objDim=2;
            prob.bounds=[0 -ones(1,prob.varDim-1);1 ones(1,prob.varDim-1)];
        case 109            % DF9
            prob.name='DF9';
            prob.func=@DF9;
            prob.varDim=10;
            prob.objDim=2;
            prob.bounds=[0 -ones(1,prob.varDim-1);1 ones(1,prob.varDim-1)];
        case 110   %DF10
            prob.name='DF10';
            prob.func=@DF10;
            prob.varDim=10;
            prob.objDim=3;
            prob.bounds=[zeros(1,2) -ones(1,prob.varDim-2);ones(1,2) ones(1,prob.varDim-2)];
        case 111   %DF11
            prob.name='DF11';
            prob.func=@DF11;
            prob.varDim=10;
            prob.objDim=3;
            prob.bounds=[zeros(1,prob.varDim);ones(1,prob.varDim)];
        case 112 % DF12
            prob.name='DF12';
            prob.func=@DF12;
            prob.varDim=10;
            prob.objDim=3;
            prob.bounds=[zeros(1,2) -ones(1,prob.varDim-2);ones(1,2) ones(1,prob.varDim-2)];
        case 113 %DF13 ->time consuming
            prob.name='DF13';
            prob.func=@DF13;
            prob.varDim=10;
            prob.objDim=3;
            prob.bounds=[zeros(1,2) -ones(1,prob.varDim-2);ones(1,2) ones(1,prob.varDim-2)];
        case 114 % DF14
            prob.name='DF14';
            prob.func=@DF14;
            prob.varDim=10;
            prob.objDim=3;
            prob.bounds=[zeros(1,2) -ones(1,prob.varDim-2);ones(1,2) ones(1,prob.varDim-2)];
    
             %% ----------------SDP test suite ----------------
        case 201
            varDim=10; objDim=3;
            bounds=[ones(1,objDim) zeros(1,varDim-objDim);
                4*ones(1,objDim) ones(1,varDim-objDim)];
            
            prob=SDP(1, varDim, objDim, bounds); % an instance of SDP1
            
        case 202
            varDim=10; objDim=3;
            bounds=[ones(1,objDim-1) -ones(1,varDim-objDim+1);
                4*ones(1,objDim-1) ones(1,varDim-objDim+1)];
            
            prob=SDP(2, varDim, objDim, bounds); % an instance of SDP2
        case 203
            varDim=10; objDim=3;
            bounds=[zeros(1,objDim-1) -ones(1,varDim-objDim+1);
                ones(1,varDim)];
            
            prob=SDP(3, varDim, objDim, bounds); % an instance of SDP3
        case 204
            varDim=10; objDim=3;
            bounds=[zeros(1,objDim-1) -ones(1,varDim-objDim+1);
                ones(1,varDim)];
            
            prob=SDP(4, varDim, objDim, bounds); % an instance of SDP4
        case 205
            varDim=10; objDim=3;
            bounds=[zeros(1,varDim); ones(1,varDim)];
            
            prob=SDP(5, varDim, objDim, bounds); % an instance of SDP5
        case 206 
            varDim=10; objDim=3;
            bounds=[zeros(1,varDim); ones(1,varDim)];
            
            prob=SDP(6, varDim, objDim, bounds); % an instance of SDP6
        case 207 
            varDim=10; objDim=3;
            bounds=[zeros(1,varDim); ones(1,varDim)];
            
            prob=SDP(7, varDim, objDim, bounds); % an instance of SDP7
        case 208
            varDim=10; objDim=3;
            bounds=[zeros(1,objDim-1) -ones(1,varDim-objDim+1);
                ones(1,varDim)];
            
            prob=SDP(8, varDim, objDim, bounds); % an instance of SDP8
        case 209
            varDim=10; objDim=3;
            bounds=[zeros(1,varDim); ones(1,varDim)];
            
            prob=SDP(9, varDim, objDim, bounds); % an instance of SDP9
        case 210
            varDim=10; objDim=3;
            bounds=[zeros(1,objDim-1) -ones(1,varDim-objDim+1);
                ones(1,varDim)];
            
            prob=SDP(10, varDim, objDim, bounds); % an instance of SDP10
        case 211
            varDim=10; objDim=3;
            bounds=[zeros(1,varDim); ones(1,varDim)];
            
            prob=SDP(11, varDim, objDim, bounds); % an instance of SDP11
        case 212
            varDim=10; objDim=3;
            nvar_a=10; nvar_b=20;
            bounds=[zeros(1,objDim-1) -ones(1,nvar_b-objDim+1);
                ones(1,nvar_b)];
            
            prob=SDP(12, varDim, objDim, bounds, nvar_a, nvar_b); % an instance of SDP12
            
        case 213 %?
            varDim=10; objDim=3;
            nobj_a=3; nobj_b=5; % min/max number of objectives
            bounds=[zeros(1,varDim); ones(1,varDim)];
            
            prob=SDP(13, varDim, objDim, bounds, nobj_a, nobj_b); % an instance of SDP13
            
        case 214
            varDim=10; objDim=3;
            bounds=[zeros(1,varDim); ones(1,varDim)];
            
            prob=SDP(14, varDim, objDim, bounds); % an instance of SDP14
        case 215 %? 8.0, 6.0, 4.0, 2.0
            varDim=10; objDim=3;
            bounds=[zeros(1,varDim); ones(1,varDim)];
            
            prob=SDP(15, varDim, objDim, bounds); % an instance of SDP15
            
        %% FDA test suite: 
        case 301 % FDA1
            prob.name='FDA1';
            prob.func=@FDA1;
            prob.varDim=10;
            prob.objDim=2;
            prob.bounds=[0 -ones(1,prob.varDim-1);1 ones(1,prob.varDim-1)];
        case 302 %  FDA2   将xiii看成空集
            prob.name='FDA2';
            prob.func=@FDA2;
            prob.varDim=10;
            prob.objDim=2;
            prob.bounds=[0 -ones(1,prob.varDim-1);1 ones(1,prob.varDim-1)];
        case 303 %  FDA3
            prob.name='FDA3';
            prob.func=@FDA3;
            prob.varDim=10;
            prob.objDim=2;
            prob.bounds=[zeros(1,2) -ones(1,prob.varDim-2);ones(1,2) ones(1,prob.varDim-2)];
        case 304 % FDA4
            prob.name='FDA4';
            prob.func=@FDA4;
            prob.varDim=10;
            prob.objDim=3;
            prob.bounds=[zeros(1,prob.varDim);ones(1,prob.varDim)];
        case 305   % FDA5
            prob.name='FDA5';
            prob.func=@FDA5;
            prob.varDim=10;
            prob.objDim=3;
            prob.bounds=[zeros(1,prob.varDim);ones(1,prob.varDim)];

        %% DMOP test suite from dCOEA
        case 401 % dmop1
            prob.name='DMOP1';
            prob.func=@DMOP1;
            prob.varDim=10;
            prob.objDim=2;
            prob.bounds=[zeros(1,prob.varDim);ones(1,prob.varDim)];
        case 402 %%%%dmop2 
            prob.name='DMOP2';
            prob.func=@DMOP2;
            prob.varDim=10;
            prob.objDim=2;
            prob.bounds=[zeros(1,prob.varDim);ones(1,prob.varDim)];
        case 403    % dMOP3
            prob.name='DMOP3';
            prob.func=@DMOP3;
            prob.varDim=10;
            prob.objDim=2;
            prob.bounds=[zeros(1,prob.varDim);ones(1,prob.varDim)];
            
       %% F test suite: A population prediction strategy for evolutionary dynamic multiobjective optimization
       % F5-F10 are defined at the end of this file, below shows an example of F8
       case 508   %%%F8
            prob.name='F8';
            prob.func=@F8;
            prob.varDim=10;
            prob.objDim=3;
            prob.bounds=[zeros(1,2) -ones(1,prob.varDim-2);ones(1,2) 2*ones(1,prob.varDim-2)];
        otherwise
            error('Problem ID is out of range or undefined. Please check problem definition')
    end

    % set the lower/upper bounds of objective/decision dimensionality for no-SDP problems
    if ~ismember(prob.name, arrayfun(@(i) ['SDP', num2str(i)], 1:15, 'UniformOutput', false))
        % the range of objective dimensionality
        prob.m_a=prob.objDim;
        prob.m_b=prob.objDim;

        % the range of decision dimensionality
        prob.n_a=prob.varDim;
        prob.n_b=prob.varDim;
    end
end


%% function to calculate time instant
function t=cal_time(tau, taut, nt, T0)
% --------------------------------------------------------|%
% The "time" term in the test suite is defined as:        |
%          t=1/nt*floor(tau/taut)                         |
% where - nt:    severity of change                       |
%       - taut:  frequency of change                      |
%       - tau:   current generation counter               |
% --------------------------------------------------------|%
% the first change occurs after T0 generations, that is, the
%generation at which a change occurs is (T0+1), (T0+taut+1), etc.  

    tau_tmp=max(tau+taut-(T0+1),0);
    t=1/nt*floor(tau_tmp/taut);
end


%% test problem definition

% Note: SDP Test Suite is defined in a separate file.

%% DF Test Suite
function xf= DF1(x,t)
    f=[];
    G=abs(sin(0.5*pi*t));
    H=0.75*sin(0.5*pi*t)+1.25;
    g=1+sum((x(:,2:end)-G).^2,2);
    f1=x(:,1);
    f2=g.*(1-(x(:,1)./g).^H);
    f=[f1,f2];
    xf=[x,f];
end

function xf=DF2(x,t)
    f=[];
    n=size(x,2);
    G=abs(sin(0.5*pi*t));
    r=1+floor((n-1)*G);
    tmp=setdiff(1:n,r);
    g=1+sum((x(:,tmp)-G).^2,2);
    f1=x(:,r);
    f2=g.*(1-(x(:,r)./g).^0.5);
    f=[f1,f2];
    xf=[x, f];
end

function xf= DF3(x,t)
    f=[];
    G=sin(0.5*pi*t);
    H=G+1.5;
    g=1+sum((x(:,2:end)-G-x(:,1).^H).^2,2);
    f1=x(:,1);
    f2=g.*(1-(x(:,1)./g).^H);
    f=[f1,f2];
    xf=[x,f];
end

function xf= DF4(x,t)
    f=[];
    n=size(x,2);

    a=sin(0.5*pi*t);
    b=1+abs(cos(0.5*pi*t));
    c=max(abs(a), a+b);
    H=1.5+a;
    g=1+sum((x(:,2:end)-a*(x(:,1)/c).^2./[2:n]).^2,2);
    f1=g.*abs(x(:,1)-a).^H;
    f2=g.*abs(x(:,1)-a-b).^H;
    f=[f1,f2];
    xf=[x,f];
end

function xf= DF5(x,t)
    f=[];
    G=sin(0.5*pi*t);
    w=floor(10*G);
    g=1+sum((x(:,2:end)-G).^2,2);
    f1=g.*(x(:,1)+0.02*sin(w*pi*x(:,1)));
    f2=g.*(1-x(:,1)+0.02*sin(w*pi*x(:,1)));
    f=[f1,f2];
    xf=[x,f];
end

function xf= DF6(x,t)
    f=[];
    G=sin(0.5*pi*t);
    a=0.2+2.8*abs(G);
    y=x(:,2:end)-G;
    g=1+sum((abs(G)*y.^2-10*cos(2*pi*y)+10),2);
    f1=g.*(x(:,1)+0.1*sin(3*pi*x(:,1))).^a;
    f2=g.*(1-x(:,1)+0.1*sin(3*pi*x(:,1))).^a;
    f=[f1,f2];
    xf=[x,f];
end

function xf= DF7(x,t)
    f=[];
    a=5*cos(0.5*pi*t);
    tmp=1./(1+exp(a*(x(:,1)-2.5)));
    g=1+sum((x(:,2:end)-tmp).^2,2);
    f1=g.*(1+t)./x(:,1);
    f2=g.*x(:,1)/(1+t);
    f=[f1,f2];
    xf=[x,f];
end

function xf= DF8(x,t)
    f=[];
    G=sin(0.5*pi*t);
    a=2.25+2*cos(2*pi*t);
    b=100*G^2;
    tmp=G*sin(4*pi*x(:,1).^b)/(1+abs(G));
    g=1+sum((x(:,2:end)-tmp).^2,2);
    f1=g.*(x(:,1)+0.1*sin(3*pi*x(:,1)));
    f2=g.*(1-x(:,1)+0.1*sin(3*pi*x(:,1))).^a;
    f=[f1,f2];
    xf=[x,f];
end

function xf= DF9(x,t)
    f=[];
    N=1+floor(10*abs(sin(0.5*pi*t)));
    g=ones(size(x,1),1);
    for i=2:size(x,2)
        tmp=x(:,i)-cos(4*t+x(:,1)+x(:,i-1));
        g=g+tmp.^2;
    end
    f1=g.*(x(:,1)+max(0, (0.1+0.5/N).*sin(2*N*pi*x(:,1))));
    f2=g.*(1-x(:,1)+max(0, (0.1+0.5/N).*sin(2*N*pi*x(:,1))));
    f=[f1,f2];
    xf=[x,f];
end

function xf= DF10(x,t)
    f=[];
    G=sin(0.5*pi*t);
    H=2.25+2*cos(0.5*pi*t);
    tmp=sin(2*pi*(x(:,1)+x(:,2)))/(1+abs(G));
    g=1+sum((x(:,3:end)-tmp).^2,2);
    f1=g.*sin(0.5*pi*x(:,1)).^H;
    f2=g.*sin(0.5*pi*x(:,2)).^H.*cos(0.5*pi*x(:,1)).^H;
    f3=g.*cos(0.5*pi*x(:,2)).^H.*cos(0.5*pi*x(:,1)).^H;
    f=[f1,f2,f3];
    xf=[x,f];
end

function xf= DF11(x,t)
    f=[];
    G=abs(sin(0.5*pi*t));
    g=1+G+sum((x(:,3:end)-0.5*G*x(:,1)).^2,2);
    y=pi*G/6+(pi/2-pi*G/3)*x(:,1:2);
    f1=g.*sin(y(:,1));
    f2=g.*sin(y(:,2)).*cos(y(:,1));
    f3=g.*cos(y(:,2)).*cos(y(:,1));
    f=[f1,f2,f3];
    xf=[x,f];
end

function xf= DF12(x,t)
    f=[];
    k=10*sin(pi*t);
    tmp1=x(:,3:end)-sin(t*x(:,1));
    tmp2=sin(floor(k*(2*x(:,1:2)-1))*pi/2);
    g=1+sum(tmp1.^2,2)+prod(tmp2,2);
    f1=g.*cos(0.5*pi*x(:,2)).*cos(0.5*pi*x(:,1));
    f2=g.*sin(0.5*pi*x(:,2)).*cos(0.5*pi*x(:,1));
    f3=g.*sin(0.5*pi*x(:,1));
    f=[f1,f2,f3];
    xf=[x,f];
end

function xf= DF13(x,t)
    f=[];
    G=sin(0.5*pi*t);
    p=floor(6*G);
    g=1+sum((x(:,3:end)-G).^2,2);
    f1=g.*cos(0.5*pi*x(:,1)).^2;
    f2=g.*cos(0.5*pi*x(:,2)).^2;
    f3=g.*sin(0.5*pi*x(:,1)).^2+sin(0.5*pi*x(:,1)).*cos(p*pi*x(:,1)).^2 +...
        sin(0.5*pi*x(:,2)).^2+sin(0.5*pi*x(:,2)).*cos(p*pi*x(:,2)).^2;
    f=[f1,f2,f3];
    xf=[x,f];
end

function xf= DF14(x,t)
    f=[];
    G=sin(0.5*pi*t);
    g=1+sum((x(:,3:end)-G).^2,2);
    y=0.5+G*(x(:,1)-0.5);
    f1=g.*(1-y+0.05*sin(6*pi*y));
    f2=g.*(1-x(:,2)+0.05*sin(6*pi*x(:,2))).*(y+0.05*sin(6*pi*y));
    f3=g.*(x(:,2)+0.05*sin(6*pi*x(:,2))).*(y+0.05*sin(6*pi*y));
    f=[f1,f2,f3];
    xf=[x,f];
end

%% FDA Test Suite
function xf= FDA1(x,t)
    f=[];
    G=sin(0.5*pi*t);
    g=1+sum((x(:, 2:end)-G).^2, 2);
    f1=x(:,1);
    f2=g.*(1-sqrt(x(:,1)./g));
    f=[f1,f2];
    xf=[x,f];
 end

function xf= FDA2(x,t)
    f=[];
    n=size(x,2);

    f1=x(:,1);
    H=2*sin(0.5*pi*(t-1));
    temp=x(:,2:n-7);
    g=1+sum(temp.^2, 2);

    ntemp=x(:,n-8:end);
    sybc=H+sum((ntemp-H/4).^2, 2);
    arbit=2^sybc;

    f2=g.*(1-(f1./g).^arbit);
    f=[f1,f2];
    xf=[x,f];
end

function xf= FDA3(x,t)
    f=[];
    n=size(x,2);

    G=abs(sin(0.5*pi*t));
    F=10^(2*sin(0.5*pi*t));

    f1=0.5*(x(:,1)^F+x(:,2)^F);

    temp=x(:,3:n);
    g=(1+G)+sum((temp-G).^2, 2);

    f2=g.*(1-(f1./g).^0.5);
    f=[f1,f2];
    xf=[x,f];
end

 function xf= FDA4(x,t)
     f=[];
     n=size(x,2);

     G=sin(0.5*pi*t);
     g=sum((x(:,3:n)-G).^2, 2);


     f1=(1+g).*cos(0.5*pi*x(:,1)).*cos(0.5*pi*x(:,2));
     f2=(1+g).*cos(0.5*pi*x(:,1)).*sin(0.5*pi*x(:,2));
     f3=(1+g).*sin(0.5*pi*x(:,1));

     f=[f1,f2,f3];
     xf=[x,f];
 end

 function xf= FDA5(x,t)
     f=[];
     G=abs(sin(0.5*pi*t));
     F=1+100*sin(0.5*pi*t)^4;

     temp=x(:,3:n);
     g=G+sum((temp-G).^2, 2);


     f1=(1+g).*cos(0.5*pi*(x(:,1).^F)).*cos(0.5*pi*(x(:,2).^F));
     f2=(1+g).*cos(0.5*pi*(x(:,1).^F)).*sin(0.5*pi*(x(:,2).^F));
     f3=(1+g).*sin(0.5*pi*(x(:,1).^F));

     f=[f1,f2,f3];
     xf=[x,f];
 end

%% DMOP Test Suite
function xf= DMOP1(x,t)
    f=[];
    n=size(x,2);
    H=0.75*sin(0.5*pi*t)+1.25;
    g=1+9*sum((x(:,2:end)).^2, 2)/(n-1);

    f1=x(:,1);
    f2=g.*(1-(x(:,1)./(g.^H)));
    f=[f1,f2];
    xf=[x,f];
end

function xf= DMOP2(x,t)
    f=[];
    G=abs(sin(0.5*pi*t));
    H=0.75*sin(0.5*pi*t)+1.25;
    g=1+sum((x(:,2:end)-G).^2, 2);
    
    f1=x(:,1);
    f2=g.*(1-(x(:,1)./g)^H);
    f=[f1,f2];
    xf=[x,f];
end

function xf= DMOP3(x,t)
    f=[];
    n=size(x,2);
    G=abs(sin(0.5*pi*t));
    
    r=randperm(n,1); % set a random variable as PF shape-related variable

    tmp=setdiff(1:n,r);
    g=1+sum((x(:,tmp)-G).^2, 2);
    
    f1=x(:,r);
    f2=g.*(1-(x(:,r)./g)^0.5);
    f=[f1,f2];
    xf=[x,f];
end

%%% F test suite: A population prediction strategy for evolutionary dynamic multiobjective optimization
%% please note: the F suite is not tested yet

function xf= F5(x,t) %%%F5
f=[];
n=size(x,2);
H=1.25+0.75*sin(pi*t);
a=2*cos(pi*t)+2;
b=2*sin(2*pi*t)+2;
indexD=1:n;
I1=indexD(find(mod(indexD,2)==0));
I2=indexD(find(mod(indexD,2)~=0));
for k=1:n
    m=H+k/n;
    qq=(abs(x(1)-a))^m;
    y(k)=x(k)-b-1+qq;
end

W=sum(y(I1).^2, 2);
f1=abs(x(1)-a)^H+W;
WW=sum(y(I2).^2, 2);
f2=abs(x(1)-a-1)^H+WW;

f=[f1,f2];
xf=[x,f];
end

function xf= F6(x,t)   %%%F6
    f=[];
    n=size(x,2);
    
    H=1.25+0.75*sin(pi*t);
    a=2*cos(1.5*pi*t)*sin(0.5*pi*t)+2;
    b=2*cos(1.5*pi*t)*cos(0.5*pi*t)+2;
    indexD=1:n;
    I1=indexD(find(mod(indexD,2)==0));
    I2=indexD(find(mod(indexD,2)~=0));
    
    for k=1:n
        m=H+k/n;
        qq=(abs(x(1)-a))^m;
        y(k)=x(k)-b-1+qq;
    end
    
    W=sum(y(I1).^2, 2);
    f1=abs(x(1)-a)^H+W;
    WW=sum(y(I2).^2, 2);
    f2=abs(x(1)-a-1)^H+WW;
    
    f=[f1,f2];
    xf=[x,f];
end

function xf= F7(x,t) %%%F7
    f=[];
    n=size(x,2);

    H=1.25+0.75*sin(pi*t);
    
    a=1.7*(1-sin(pi*t))*sin(pi*t)+3.4;
    b=1.4*(1-sin(pi*t))*cos(pi*t)+2.1;
    
    indexD=1:n;
    I1=indexD(find(mod(indexD,2)==0));
    I2=indexD(find(mod(indexD,2)~=0));
    
    for k=1:n
        m=H+k/n;
        qq=(abs(x(1)-a))^m;
        y(k)=x(k)-b-1+qq;
    end
    W=sum(y(I1).^2, 2);
    f1=abs(x(1)-a)^H+W;
    WW=sum(y(I2).^2, 2);
    f2=abs(x(1)-a-1)^H+WW;
    
    f=[f1,f2];
    xf=[x,f];
end

function xf= F8(x,t)     %%%F8
     f=[];
     G=sin(0.5*pi*t);
     H=1.25+0.75*sin(pi*t);
     gg=((x(:,1)+x(:,2))/2).^H;
     g=sum((x(:,3:end)-gg-G).^2,2);

     f1=(1+g).*cos(0.5*pi*x(:,2)).*cos(0.5*pi*x(:,1));
     f2=(1+g).*cos(0.5*pi*x(:,2)).*sin(0.5*pi*x(:,1));
     f3=(1+g).*sin(0.5*pi*x(:,2));

     f=[f1,f2,f3];
     xf=[x,f];
end
 
function xf= F9(x,t)   %%%F9
    f=[];
    n=size(x,2);
    
    H=1.25+0.75*sin(pi*t);
    
    a=2*cos((t-floor(t))*pi)+2;
    b=2*sin(2*(t-floor(t))*pi)+2;
    
    indexD=1:n;
    I1=indexD(find(mod(indexD,2)==0));
    I2=indexD(find(mod(indexD,2)~=0));
    
    for k=1:n
        m=H+k/n;
        qq=(abs(x(1)-a))^m;
        y(k)=x(k)-b-1+qq;
    end
    W=sum(y(I1).^2, 2);
    f1=abs(x(1)-a)^H+W;
    WW=sum(y(I2).^2, 2);
    f2=abs(x(1)-a-1)^H+WW;
    
    f=[f1,f2];
    xf=[x,f];
end

function xf= F10(x,t)   %%%F10
    f=[];
    n=size(x,2);
    
    H=1.25+0.75*sin(pi*t);
    
    a=2*cos(pi*t)+2;
    b=2*sin(2*pi*t)+2;
    
    indexD=1:n;
    I1=indexD(find(mod(indexD,2)==0));
    I2=indexD(find(mod(indexD,2)~=0));
    
    for k=1:n
        m=H+k/n;
        qq=(abs(x(1)-a))^m;
        if mod(tau,2)~=0 
            y(k)=x(k)-b-qq;
        else
            y(k)=x(k)-b-1+qq;
        end
    end
    W=sum(y(I1).^2, 2);
    f1=abs(x(1)-a)^H+W;
    WW=sum(y(I2).^2, 2);
    f2=abs(x(1)-a-1)^H+WW;
    
    f=[f1,f2];
    xf=[x,f];
end