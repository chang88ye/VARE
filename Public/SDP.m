classdef SDP<handle
    % Base class of SDP test suite

	properties(SetAccess = protected)
        name;           % probblem name
        func;           % problem handle
        
		varDim;         % the actual number of variables
		objDim;         % the actual number of objectives
        
        n_a;            % the lower bound of varDim
        n_b;            % the upper bound of varDim
        
        bounds;         % lower and upper boudns of variables
        
        m_a;            % the lower bound of objDim
        m_b;            % the upper bound of objDim

        v_;             % current varying parameters in SDP problems
        u_;             % storing the previous v_

        curT;           % current time point
        
	end

	methods
        %% Constructor

		function prob=SDP(id, nvar, mobj, bounds, varargin)
            % varagin stores either bounds of objective dimensions or bounds of variable dimensions
            prob.varDim=nvar;
            prob.objDim=mobj;
            prob.bounds=bounds;
            prob.v_=[];
            prob.name =['SDP', num2str(id)];
            
            prob.func =@(x,t, varargin)prob.(prob.name)(x,t, varargin); %str2func (['SDP', num2str(id)]);
            
            % assign intial values to dynamic parameters stored in prob.v_
            switch id
                case 1
                    prob.v_=(mobj+1:nvar)/nvar;
                case 4
                    prob.v_=1;
                case 7
                    prob.v_=[0.5, 2:5];
                case 12
                    assert(nargin==6);
                    assert(nvar>=varargin{1} && nvar<=varargin{2})
                    
                    prob.n_a=varargin{1};
                    prob.n_b=varargin{2};
                    prob.v_=nvar;
                case 13
                    assert(nargin==6);
                    assert(mobj>=varargin{1} && mobj<=varargin{2})
                    
                    prob.m_a=varargin{1};
                    prob.m_b=varargin{2};
                    prob.v_=mobj;
                case 15
                    prob.v_=(mobj-1)*ones(1,2);
                otherwise
                    prob.v_=[];

                    prob.m_a=prob.objDim;
                    prob.m_b=prob.objDim;
                    
                    prob.n_a=prob.varDim;
                    prob.n_b=prob.varDim;
            end
            
            prob.curT=0;
        end
        
        function t=update_time(obj, tau, taut, nt, T0) % same as cal_time 
            tau_tmp=max(tau+taut-(T0+1),0);
            t=1/nt*floor(tau_tmp/taut);
        end
        
        function obj = updateSDP1(obj, t)
            % update y_i,t
            obj.curT=t;
            
            x_ = obj.v_ + 5* (rand(1,length(obj.v_))-0.5)*sin(0.5*pi*t);
            x_(x_>1 | x_<0)= rand(1,sum(x_>1 | x_<0));
            obj.v_ =x_;
            
        end
        
        function obj = updateSDP4(obj, t)
            % update rnd_t
            obj.curT=t;
            
            obj.v_=rand;
            
        end
        
        function obj = updateSDP7(obj, t)
            % update h
            obj.curT=t;
            
            obj.v_=1:5;
            pt=randi(5);
            obj.v_(pt)=0.5;
        end
        
        function obj = updateSDP12(obj, t)
            % update n_t
            obj.curT=t;
            
            obj.v_= randi([obj.n_a,obj.n_b],1);
            obj.varDim=obj.v_;
        end
        
        function obj = updateSDP13(obj, t)
            % update m
            obj.curT=t;
            
            obj.v_= randi([obj.m_a ,obj.m_b],1); 
            obj.objDim=obj.v_;
        end
        
        function obj = updateSDP15(obj, t)
            % update dt, pt
            
            obj.curT=t;
            
            obj.v_=randi(obj.objDim-1,1,2);%[dt, pt]
        end
        
        % ======function evaluation======= %
        function xf=SDP1(obj, x, t, varargin)
            
            if t>obj.curT % time changes
                obj.u_=obj.v_;
                obj=updateSDP1(obj, t);
            end
            
            if ~isempty(varargin{1})
                g=ones(size(x,1),1);
            else
                g=1+sum((x(:,obj.objDim+1:end)-obj.v_).^2,2);
            end
            
            f=g.*(x(:,1:obj.objDim)./repmat(prod(x(:,obj.objDim),2), 1, obj.objDim)).^0.5;

            if ~isempty(varargin{1})
                xf=f;
            else
                xf=[x, f];
            end
        end
        
        function xf=SDP2(obj, x, t, varargin)
            if t>obj.curT % time changes
                obj.u_=obj.v_;
                obj.curT=t;
            end
            
            nM=obj.objDim-1;
            
            tmp=cos(t+2*x(:,1));
            
            if ~isempty(varargin{1})
                g=ones(size(x,1),1);
            else
                g=1+sin(pi*x(:,1)/8).*sum((x(:,nM+1:end)-tmp).^2, 2);
            end
            
            fi=(1+t+sum(x(:,1:nM),2))./x(:,1:nM)-1.0;
            fM=sum(x(:,1:nM),2)/(1+t);
            f=g.*[fi, fM];       

            if ~isempty(varargin{1})
                xf=f;
            else
                xf=[x, f];
            end
            
        end
        
        function xf=SDP3(obj, x, t, varargin)
            if t>obj.curT % time changes
                obj.u_=obj.v_;
                obj.curT=t;
            end
            
            nM=obj.objDim-1;

            pt=floor(5*abs(sin(pi*t)));
            y=x(:,nM+1:end)-cos(t);
            
            if ~isempty(varargin{1})
                g=ones(size(x,1),1);
            else
                g=1+sum(4*y.^2-cos(2*pt*pi*y)+1, 2);
            end

            f=g.*cumprod([ones(size(g,1),1),x(:,1:nM)+0.05*sin(6*pi*x(:,1:nM))],2)...
                .*[1-x(:,1:nM)+0.05*sin(6*pi*x(:,1:nM)), ones(size(g,1),1)];

            if ~isempty(varargin{1})
                xf=f;
            else
                xf=[x, f];
            end
        end
        
        
        function xf=SDP4(obj, x, t, varargin)
            if t>obj.curT % time changes
                obj.u_=obj.v_;
                obj=updateSDP4(obj, t);
            end
            
            nM=obj.objDim-1;
    
            w=sign(obj.v_-0.5)*floor(6*abs(sin(0.5*pi*t)));
            s=mean(x(:, 1:nM),2);

            if ~isempty(varargin{1})
                g=ones(size(x,1),1);
            else
                g=1+sum((x(:,nM+1:end)-cos(t+x(:,1)+x(:,nM:end-1))).^2, 2);
            end
            
            fend=[s+0.05*sin(w*pi*s), 1-s+0.05*sin(w*pi*s)];
            f=g.*[x(:, 1:nM-1), fend];

            if ~isempty(varargin{1})
                xf=f;
            else
                xf=[x, f];
            end
        end
        
        
        function xf=SDP5(obj, x, t, varargin)
            if t>obj.curT % time changes
                obj.u_=obj.v_;
                obj.curT=t;
            end
            
            nM=obj.objDim-1;
    
            Gt=abs(sin(0.5*pi*t));
            y=pi*Gt+(pi/2-pi/3*Gt)*x(:,1:nM);

            if ~isempty(varargin{1})
                g=ones(size(x,1),1);
            else
                g=1+Gt+sum((x(:, nM+1:end)-0.5*Gt*x(:,1)).^2,2);
            end
            
            f=g.*cumprod([ones(size(g,1),1),cos(y(:,1:nM))],2)...
                .*[sin(y(:,1:nM)), ones(size(g,1),1)];

            if ~isempty(varargin{1})
                xf=f;
            else
                xf=[x, f];
            end
        end
        
        function xf=SDP6(obj, x,t, varargin)
            if t>obj.curT % time changes
                obj.u_=obj.v_;
                obj.curT=t;
            end
            
            nM=obj.objDim-1;
            
            a=0.5*abs(sin(pi*t));
            k=10*cos(2.5*pi*t);
            
            if ~isempty(varargin{1})
                g=ones(size(x,1),1);
            else
                g=1+sum((x(:,nM+1:end)-0.5).^2.*(1+abs(cos(8*pi*x(:,nM+1:end)))),2);
            end
            
            f=g.*fliplr(cumprod([ones(size(g,1),1),cos(0.5*pi*x(:,1:nM))],2))...
                .*[ones(size(g,1),1),sin(0.5*pi*x(:,nM:-1:1))];
            
            tmp=x(:,1)<a;
            f(tmp,end)=g(tmp).*abs(k*(cos(0.5*pi*x(tmp,1))-cos(0.5*pi*a))+sin(0.5*pi*a));
            
            if ~isempty(varargin{1})
                xf=f;
            else
                xf=[x, f];
            end
        end
        
        function xf=SDP7(obj, x,t, varargin)
            if t>obj.curT % time changes
                obj.u_=obj.v_;
                obj=obj.updateSDP7(t);
            end
            
            nM=obj.objDim-1;
            
            h=obj.v_; % valleys
            y=0:2:8;

            if ~isempty(varargin{1})
                g=ones(size(x,1),1);
            else
                hl=@(v) min(h+10*(10*v-y).^2, [], 2); % for each column vector
                g=0.5+mean(cell2mat(arrayfun(hl, x(:,nM+1:end), 'UniformOutput', false)),2);
            end
            
            f=g.*cumprod([ones(size(g,1),1),x(:,1:nM)],2)...
                .*[1-x(:,1:nM),ones(size(g,1),1)];
            
            if ~isempty(varargin{1})
                xf=f;
            else
                xf=[x, f];
            end
        end
        
        function xf=SDP8(obj,x,t, varargin)
            if t>obj.curT % time changes
                obj.u_=obj.v_;
                obj.curT=t;
            end
            
            nM=obj.objDim-1;
            
            k=10*sin(pi*t);
            
            if ~isempty(varargin{1})
                g=ones(size(x,1),1);
            else
                g=1+sum((x(:,nM+1:end)-sin(t*x(:,1))).^2,2)...
                    +abs(prod(sin(floor(k*(2*x(:,1:nM)-1))*pi/2),2));
            end
            
            f=g.*fliplr(cumprod([ones(size(g,1),1),cos(0.5*pi*x(:,1:nM))],2))...
                .*[ones(size(g,1),1),sin(0.5*pi*x(:,nM:-1:1))];
            
            if ~isempty(varargin{1})
                xf=f;
            else
                xf=[x, f];
            end
        end
        
        function xf=SDP9(obj,x,t, varargin)
            if t>obj.curT % time changes
                obj.u_=obj.v_;
                obj.curT=t;
            end
            
            nM=obj.objDim-1;
            
            G=abs(sin(0.5*pi*t));
            p=floor(6*G);
            
            if ~isempty(varargin{1})
                g=ones(size(x,1),1);
            else
                g=1+sum((x(:,nM+1:end)-abs(atan(cot(3*pi*t^2)))/pi).^2,2);
            end
            
            fm=sum(sin(0.5*pi*x(:,1:nM)).^2+sin(0.5*pi*x(:,1:nM)).*cos(p*pi*x(:,1:nM)).^2,2);
            
            f=G+[g.*cos(0.5*pi*x(:,1:nM)).^2, fm];
            
            if ~isempty(varargin{1})
                xf=f;
            else
                xf=[x, f];
            end
        end
        
        function xf=SDP10(obj, x,t, varargin)
            if t>obj.curT % time changes
                obj.u_=obj.v_;
                obj.curT=t;
            end
            
            nM=obj.objDim-1;
            
            s=mean(x(:,1:nM).^2,2);
            r=floor(10*abs(sin(0.5*pi*t)));
            
            if ~isempty(varargin{1})
                g=ones(size(x,1),1);
            else
                g=1+sum((x(:,nM+1:end)-sin(x(:,1)+0.5*pi*t)).^2, 2);
            end
            
            f=[x(:,1:nM),g.*(2-s-s.^0.5.*(-sin(2.5*pi*s)).^r)];
            
            if ~isempty(varargin{1})
                xf=f;
            else
                xf=[x, f];
            end
        end
        
        function xf=SDP11(obj, x,t, varargin)
            if t>obj.curT % time changes
                obj.u_=obj.v_;
                obj.curT=t;
            end
            
            nM=obj.objDim-1;
            
            p=x(:,nM+1:end)-abs(sin(0.5*pi*t));
            phi=sum(x(:,1:nM),2)<mod(3*t+0.2,1) & sum(x(:,1:nM),2)> mod(3*t,1);
            
            if ~isempty(varargin{1})
                g=ones(size(x,1),1);
            else
                g=phi.*sum(-0.9*p.^2+abs(p).^(0.6),2)+(1-phi).*sum(p.^2,2);
            end
            
            f=g.*cumprod([ones(size(g,1),1),cos(0.5*pi*x(:,1:nM))],2)...
                .*[sin(0.5*pi*x(:,1:nM)), ones(size(g,1),1)];
            
            if ~isempty(varargin{1})
                xf=f;
            else
                xf=[x, f];
            end
        end
        
        
        function xf=SDP12(obj, x,t, varargin)
            if t>obj.curT % time changes
                obj.u_=obj.v_;
                obj=obj.updateSDP12(t);
            end
            
            nM=obj.objDim-1;
            
            nD =obj.varDim;
            
            if ~isempty(varargin{1}) % generate PF.
                g=ones(size(x,1),1);
            else %  fitness evaluation
                %                 N=size(x,1);
                %                 if obj.u_<nD % increasing variables
                %                     pop=init_Pop(nD,N,obj.bounds); %new variables randomly initialised
                %                     x =[x, pop(:,1:nD-obj.u_)]; %new variables appended to x
                %                 else
                %                     x =x(:,1:nD); % backward remove abandunt variables
                %                     %%%% OR remove randomly some distance-related variables as follows
                %                     % xp=nM:obj.u_;
                %                     % idx= randperm(obj.u_ -nM, nD-nM)
                %                     % x=x(:, [1:nM, xp(idx)])
                %                 end

                g=1+sum((x(:,nM+1:nD)-sin(nD*t)*sin(2*pi*x(:,1))), 2);
            end
            
            f=g.*cumprod([ones(size(g,1),1),x(:,1:nM)],2)...
                .*[1-x(:,1:nM),ones(size(g,1),1)];
            
            if ~isempty(varargin{1})
                xf=f;
            else
                xf=[x, f];
            end
        end
        
        
        function xf=SDP13(obj, x,t, varargin)
            if t>obj.curT % time changes
                obj.u_=obj.v_;
                obj=obj.updateSDP13(t);
            end
            
            

            if ~isempty(varargin{1})
                nM=obj.u_ -1; % previous objective demensionality
                g=ones(size(x,1),1);

            else
                nM=obj.objDim-1;
                idx=nM+1:obj.varDim;
                g=1+sum((x(:,idx)-idx*t./(nM+idx*t)).^2,2);
            end

            y=pi*(1+x(:,1:nM))/6;
            
            f=g.*cumprod([ones(size(g,1),1),cos(y(:,1:nM))],2)...
                .*[sin(y(:,1:nM)), 0.5*ones(size(y,1),1)];
            
            if ~isempty(varargin{1})
                xf=f;
            else
                xf=[x, f];
            end
        end
        
        function xf=SDP14(obj,x,t, varargin)
            if t>obj.curT % time changes
                obj.u_=obj.v_;
                obj.curT=t;
            end
            
            nM=obj.objDim-1;
            
            d=1+floor((nM-1)*abs(sin(0.5*pi*t)));
            
            if ~isempty(varargin{1})
                g=ones(size(x,1),1);
            else
                g=1+sum((x(:,nM+1:end)-0.5).^2,2);
            end
            
            y=[x(:, 1:d), 0.5+x(:,d+1:end).*(g-1)*abs(0.5*pi*t)];
            
            f=g.*cumprod([ones(size(g,1),1),y(:,1:nM)],2)...
                .*[g-y(:,1:nM), ones(size(g,1),1)];
            
            if ~isempty(varargin{1})
                xf=f;
            else
                xf=[x, f];
            end
        end
        
        function xf=SDP15(obj,x,t, varargin)
            if t>obj.curT % time changes
                obj.u_=obj.v_;
                obj=obj.updateSDP15(t);
            end
            nM=obj.objDim-1;
            
            function f=computeObjectives(nM,v,x,t,pf)
                
                [d, p]=deal(v(1),v(2));
                
                if pf ==1 % to calculate pf
                    g=zeros(size(x,1),1);
                else
                    g=sum((x(:, nM+1:end)-d/nM).^2,2);
                end
                
                k=mod(p+(1:nM)-1, nM)+1;
                xk=x(:,k);
                y=[0.5*pi*xk(:,1:d), acos(1./(2^0.5*(1+xk(:,d:end).*abs(sin(0.5*pi*t)*g))))];
                
                g=g+1;
                
                f=g.^(1:nM+1).*cumprod([ones(size(g,1),1),cos(y(:,1:nM))],2)...
                    .*[sin(y(:,1:nM)), ones(size(g,1),1)];
            end
            if ~isempty(varargin{1})
                xf=computeObjectives(nM,obj.u_,x,t,1); % calculating pf
            else
                f=computeObjectives(nM,obj.v_,x,t,0); % calcuating f
                xf =[x, f];
            end
            
        end
    end
    
end