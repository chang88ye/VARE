function value=checkChange(func, pop, ObjNum, t, eta)

    value=false;
    
    [ps, dim]=size(pop);
    
    num_sensors= floor(eta*ps);
    
    for i=randperm(ps, num_sensors)
    
        old_obj=pop(i,dim-ObjNum+1:end);
        old_dec=pop(i,1:dim-ObjNum);

        
    
        new_x_obj=func(old_dec,t);
        new_obj=new_x_obj(dim-ObjNum+1:end);

        % either variable dimensionality or objective value changes
        if length(new_obj)~=length(old_obj) || any(abs(new_obj-old_obj)>1e-5)
            value=true;
            break
        end
    
    end
end