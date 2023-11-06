function [keep,selected]=filter_boundary(PopObjs)
    [N, ObjNum]=size(PopObjs);
    % remove boundary solutions 
    selected=[];
    z_min=min(PopObjs,[],1);
    z_objs=PopObjs-z_min;
    
    for obj=1:ObjNum
        w=1e-3*ones(1, ObjNum); w(obj)=1;
        % identify extreeme point using ASF along obj axis
        [~, idx_xtr]=min(max(z_objs./w, [], 2));
        selected = [selected,idx_xtr]; 
    end
    keep=true(N,1);
    keep(selected)=false;
    % tmp=PopObjs<min(PopObjs(selected,:),[],1);
    % keep(any(tmp,2))=false;
end