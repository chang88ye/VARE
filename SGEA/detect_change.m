function value=detect_change(ind, ins, t)
dim=size(ind,2);
ObjNum=ins.objDim;

value=false;
    
old_obj=ind(1,dim-ObjNum+1:end);
old_dec=ind(1,1:dim-ObjNum);

new_obj=ins.func(old_dec,t);
if any(abs(new_obj(dim-ObjNum+1:end)-old_obj)>1e-5)
    value=true;
end

end