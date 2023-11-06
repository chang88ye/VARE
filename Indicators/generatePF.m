function PF=generatePF(prob,t)

no=10000;  % expected number of PF points

switch prob.name
    case 'DF1' %101    %%%DF1
        x1=0:1/no:1;
        H=0.75*sin(0.5*pi*t)+1.25;
        f1=x1;
        f2=1-x1.^H;
        PF=[f1;f2]';
    case 'DF2' %102    %%%DF2
        x1=0:1/no:1;
        f1=x1;
        f2=1-x1.^0.5;
        PF=[f1;f2]';
    case 'DF3' %  %%%DF3
        x1=0:1/no:1;
        H=1.5+sin(0.5*pi*t);
        f1=x1;
        f2=1-x1.^H;
        PF=[f1;f2]';
    case 'DF4' %   %%%DF4
        bt=1+abs(cos(0.5*pi*t));
        H=1.5+sin(0.5*pi*t);
        x1=0:1/no:bt;
        f1=x1.^H;
        f2=(bt-x1).^H;
        PF=[f1;f2]';
    case 'DF5' %  %%%DF5
        x1=0:1/no:1;
        w=floor(10*sin(0.5*pi*t));
        f1=(x1+0.02*sin(w*pi*x1));
        f2=(1-x1+0.02*sin(w*pi*x1));
        PF=[f1;f2]';
    case 'DF6' %  %%%DF6
        x1=0:1/no:1;
        at=0.2+2.8*abs(sin(0.5*pi*t));
        f1=(x1+0.1*sin(3*pi*x1)).^at;
        f2=(1-x1+0.1*sin(3*pi*x1)).^at;
        PF=[f1;f2]';
    case 'DF7' %  %%%DF7
        x1=1:1/no:4;
        at=1+t;
        f1=at./x1;
        f2=x1./at;
        PF=[f1;f2]';
    case 'DF8' %  %%%DF8
        x1=0:1/no:1;
        at=2.25+2*cos(2*pi*t);
        f1=(x1+0.1*sin(3*pi*x1));
        f2=(1-x1+0.1*sin(3*pi*x1)).^at;
        PF=[f1;f2]';
    case 'DF9' %  %%%DF9 - disconnected
        x1=(0:1/no:1);
        nt=1+floor(10*abs(sin(0.5*pi*t)));
        f1=x1+max(0, (0.1+0.5/nt)*sin(2*nt*pi*x1));
        f2=1-x1+max(0, (0.1+0.5/nt)*sin(2*nt*pi*x1));
        k=f1+f2>1;
        f1(k)=[];
        f2(k)=[];
        PF=[f1;f2]';
    case 'DF10' % %%%DF10
        [x1,x2]=meshgrid(linspace(0,1,100));
        h=2.25+2*cos(0.5*pi*t);
        f1=sin(0.5*pi*x1).^h;
        f2=sin(0.5*pi*x2).^h.*cos(0.5*pi*x1).^h;
        f3=cos(0.5*pi*x2).^h.*cos(0.5*pi*x1).^h;
        PF=[f1(:),f2(:),f3(:)]; % column-wise concatenation by (:)
    case 'DF11' %11 %%%DF11
        [x1,x2]=meshgrid(linspace(0,1,100));
        gt=abs(sin(0.5*pi*t));
        y1=pi/6*gt+(pi/2-pi/3*gt)*x1;
        y2=pi/6*gt+(pi/2-pi/3*gt)*x2;
        f1=sin(y1);
        f2=sin(y2).*cos(y1);
        f3=cos(y2).*cos(y1);
        PF=[f1(:),f2(:),f3(:)];
    case 'DF12' %12 %%%DF12 -disconnected

        [x1,x2]=meshgrid(linspace(0,1,100));
        k=10*sin(pi*t);
        tmp2=abs(sin(floor(k*(2*x1-1))*pi/2).*sin(floor(k*(2*x2-1))*pi/2));
        g=1+tmp2;
        f1=g.*cos(0.5*pi*x2).*cos(0.5*pi*x1);
        f2=g.*sin(0.5*pi*x2).*cos(0.5*pi*x1);
        f3=g.*sin(0.5*pi*x1);

        P=[f1(:),f2(:),f3(:)];
        [PF,~]=NondomSort(P,0,3,3);
    case 'DF13' %13 %%%DF13 -disconnected
        [x1, x2]=meshgrid(linspace(0,1,100));
        G=sin(0.5*pi*t);
        p=floor(6*G);
        f1=cos(0.5*pi*x1).^2;
        f2=cos(0.5*pi*x2).^2;
        f3=(sin(0.5*pi*x1).^2+sin(0.5*pi*x1).*cos(p*pi*x1).^2 +...
            sin(0.5*pi*x2).^2+sin(0.5*pi*x2).*cos(p*pi*x2).^2);
        P=[f1(:),f2(:),f3(:)];
        [PF,~]=NondomSort(P,0,3,3);
    case 'DF14' % %%%DF14 -degenerate
        [x1,x2]=meshgrid(linspace(0,1,100));
        y1=0.5+sin(0.5*pi*t)*(x1-0.5);
        f1=(1-y1+0.05*sin(6*pi*y1));
        f2=(1-x2+0.05*sin(6*pi*x2)).*(y1+0.05*sin(6*pi*y1));
        f3=(x2+0.05*sin(6*pi*x2)).*(y1+0.05*sin(6*pi*y1));
        PF=[f1(:),f2(:),f3(:)];
        
        %% SDP Test Suite: its evaluation functions has an option to generate PF from sampling of PF-related variables
        %
    case arrayfun(@(i) ['SDP', num2str(i)], 1:15, 'UniformOutput', false) % equivalent to {'SDP1', ..., 'SDP15'}
        if strcmp(prob.name, 'SDP1')
            nd=prob.objDim;
        else
            nd=prob.objDim-1;
        end

        if nd<2 % low dimensional case
            no=1000;
        else
            no=max(5,ceil(10000^(1/nd))); % creating around 1e4 samples regardless of objective dimensionality
        end

        if strcmp(prob.name, 'SDP1')|| strcmp(prob.name, 'SDP2')
            range=linspace(1,4,no); % PF-related variables in the range [1,4]
        else
            range=linspace(0,1,no); % PF-related variables in the range [0,1]
        end
        
        out = cell(1,nd); % create a cell array to store a comma-separated list with variable number of items

        [out{:}]=ndgrid(range); % nd output variables
        x=cellfun(@(e) e(:), out, 'UniformOutput', false);% a cell array of column vectors %vertcat(out{:});
        PF=prob.func(horzcat(x{:}),t, 1); % convert cell array to a matrix before passing into the evaluation function
        
        if strcmp(prob.name, 'SDP9')|| strcmp(prob.name, 'SDP10')|| strcmp(prob.name, 'SDP11') % disconnected problems
            PF=NondomSort(PF,0,prob.objDim,prob.objDim); % get nondominated set
        end
        %
        %% SDP Test Suite: examples of generating PF for three objectives only
        %     case 'SDP1' % 201
        %         [x,y,z]=ndgrid(1:0.2:4);
        %         f1=x./(y.*z).^0.5;
        %         f2=y./(x.*z).^0.5;
        %         f3=z./(x.*y).^0.5;
        %                
        %         PF=[f1(:),f2(:),f3(:)];
        %     case 'SDP2' %202
        %         [x,y]=meshgrid(1:0.1:4);
        %         f1=(1+y)./x;
        %         f2=(1+x)./y;
        %         f3=(x+y);
        %         
        %         PF=[f1(:),f2(:),f3(:)];
        %     case 'SDP3' %203
        %         [x,y]=meshgrid(0:0.02:1);
        %         w=6;A=0.05;
        %         f1=1-x+A*sin(w*pi*x);
        %         f2=(1-y+A*sin(w*pi*y)).*(x+A*sin(w*pi*x));
        %         f3=(x+A*sin(w*pi*x)).*(y+A*sin(w*pi*y));
        %         
        %         PF=[f1(:),f2(:),f3(:)];
        %     case 'SDP4' %204
        %         [x,y]=meshgrid(0:0.02:1);
        %         w=6;A=0.05;
        %         f1=x;
        %         s=(x+y)/2;
        %         f2=(s+A*sin(w*pi*s));
        %         f3=(1-s+A*sin(w*pi*s));
        %         
        %         PF=[f1(:),f2(:),f3(:)];
        %     case 'SDP5' %205
        %         [x1,y1]=meshgrid(0:0.02:1);
        %         g=abs(sin(0.5*pi*t));
        %         x=pi*g/6+(0.5*pi-pi*g/3)*x1;
        %         y=pi*g/6+(0.5*pi-pi*g/3)*y1;
        %         f1=(1+g).*sin(x);
        %         f2=(1+g).*sin(y).*cos(x);
        %         f3=(1+g).*cos(y).*cos(x);
        %         
        %         PF=[f1(:),f2(:),f3(:)];
        %     case 'SDP6' %206
        %         [x1,y1]=meshgrid(0:0.02:1);
        %         k=0.5*abs(sin(pi*t));
        %         p=10*cos(2.5*pi*t);
        %         x=(0.5*pi)*x1;
        %         y=(0.5*pi)*y1;
        %         f1=cos(y).*cos(x);
        %         f2=sin(y).*cos(x);
        % 
        %         tmp=(x1<k).*(sin(pi*k/2)+p*(cos((0.5*pi)*x1)-cos(pi*k/2)));
        %         tmp(tmp<0)=nan;
        %         tmp=tmp*(p<=0);
        %         f3=sin((0.5*pi)*x1).*(x1>=k)+tmp;
        %         
        %         PF=[f1(:),f2(:),f3(:)];
        %       
        %     case 'SDP10' %210
        %         [x,y]=meshgrid(0:0.2:1);
        %         
        %         s=(x.^2+y.^2)/2;
        %         p=floor(10*abs(sin(0.5*pi*t)));
        %         
        %         f1=x;
        %         f2=y;
        %         f3=2-(s+s.^0.5.*(-1*sin(2.5*pi*s)).^p);
        %         
        %         PF=eff_pf(f1,f2,f3);
        % 
        %     case 'SDP11' %211
        %         x=0:0.001:1;
        %         f1=(1-x);
        %         f2=0.5*x;
        %         f3=0.5*x;
        %         
        %         PF=[f1(:),f2(:),f3(:)];
        %     case 'SDP12' %212
        %         x=0:0.001:1;
        %         f1=sin(0.5*pi*x);
        %         f2=cos(0.5*pi*x)/2.^0.5;
        %         f3=cos(0.5*pi*x)/2.^0.5;
        % 
        %         PF=[f1(:),f2(:),f3(:)];
        %     case 'SDP13' %213
        %         [x1,y1]=meshgrid(0:0.02:1);
        %         x=pi*(1+x1)/6;
        %         y=pi*(1+y1)/6;
        %         f1=sin(x);
        %         f2=sin(y).*cos(x);
        %         f3=sin(pi*(1+2/3)/6).*cos(y).*cos(x);
        % 
        %         PF=[f1(:),f2(:),f3(:)];
    
        
        %% FDA Test Suite
    case 'FDA1' %301   %%%fad1
        x1=0:1/no:1;
        f1=x1;
        f2=1-x1.^0.5;
        PF=[f1;f2]';
    case 'FDA2' %302  %%%fad2
        x1=0:1/no:1;
        h=1.5+sin(0.5*pi*t);
        H=2*h;
        f1=x1;
        f2=1-x1.^H;
        PF=[f1;f2]';
    case 'FDA3' %303  %%%fad3
        x1=0:1/no:1;
        x2=0:1/no:1;
        gt=abs(sin(0.5*pi*t));
        ft=10^(2*sin(0.5*pi*t));
        f1=0.5*(x1.^ft+x2.^ft);
        f2=(1+gt).*(1-sqrt(f1./(1+gt)));
        PF=[f1;f2]';
    case 'FDA4' %304 %%%FAD4
        [x1,x2]=meshgrid(linspace(0,1,100));
        f1=cos(0.5*pi*x1).*cos(0.5*pi*x2);
        f2=cos(0.5*pi*x1).*sin(0.5*pi*x2);
        f3=sin(0.5*pi*x1);
        PF=[f1(:),f2(:),f3(:)];
    case 'FDA5' %305 %%%FAD5
        [x1,x2]=meshgrid(linspace(0,1,100));
        gt=abs(sin(0.5*pi*t));
        ft=1+100*(sin(0.5*pi*t)^4);
        f1=(1+gt).*cos(0.5*pi*x1.^ft).*cos(0.5*pi*x2.^ft);
        f2=(1+gt).*cos(0.5*pi*x1.^ft).*sin(0.5*pi*x2.^ft);
        f3=(1+gt).*sin(0.5*pi*x1.^ft);
        PF=[f1(:),f2(:),f3(:)];
        
        %% DMOP Test Suite
    case 'DMOP1' %401   %%%dmop1
        x1=0:1/no:1;
        H=0.75*sin(0.5*pi*t)+1.25;
        f1=x1;
        f2=1-x1.^H;
        PF=[f1;f2]';
    case 'DMOP2' %402   %%%dmop2
        x1=0:1/no:1;
        H=0.75*sin(0.5*pi*t)+1.25;
        f1=x1;
        f2=1-x1.^H;
        PF=[f1;f2]';
    case 'DMOP3' %403   %%%dmop3
        x1=0:1/no:1;
        f1=x1;
        f2=1-x1.^0.5;
        PF=[f1;f2]';
        
        %% F Test Suite
    case 'F5' %505   %%%F5
        x1=0:1/no:1;
        H=1.25+sin(0.5*pi*t)*0.75;
        f1=x1.^H;
        f2=(1-x1).^H;
        PF=[f1;f2]';
    case 'F6' %506   %%%F6
        x1=0:1/no:1;
        H=1.25+sin(0.5*pi*t)*0.75;
        f1=x1.^H;
        f2=(1-x1).^H;
        PF=[f1;f2]';
    case 'F7' %507   %%%F7
        x1=0:1/no:1;
        H=1.25+sin(0.5*pi*t)*0.75;
        f1=x1.^H;
        f2=(1-x1).^H;
        PF=[f1;f2]';
    case 'F8' %508 %%%F8
        [x1,x2]=meshgrid(linspace(0,1,100));
        f1=cos(0.5*pi*x1).*cos(0.5*pi*x2);
        f2=cos(0.5*pi*x1).*sin(0.5*pi*x2);
        f3=sin(0.5*pi*x1);
        PF=[f1(:),f2(:),f3(:)];
    case 'F9' %509   %%%F9
        x1=0:1/no:1;
        H=1.25+sin(0.5*pi*t)*0.75;
        f1=x1.^H;
        f2=(1-x1).^H;
        PF=[f1;f2]';
    case 'F10' %510   %%%F10
        x1=0:1/no:1;
        H=1.25+sin(0.5*pi*t)*0.75;
        f1=x1.^H;
        f2=(1-x1).^H;
        PF=[f1;f2]';
end
end



%% helper functions

function f=eff_pf(f1,f2,f3) % get efficient front

h =[f1(:),f2(:),f3(:)];

f=NondomSort(h,0,size(h,2), size(h,2));

end


