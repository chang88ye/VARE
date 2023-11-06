function newx=checkbound(newx,lu)
[m,n]=size(newx);
for i=1:m
    for j=1:n
        if newx(i,j)>lu(2,j)
            newx(i,j)=2*lu(2,j)-newx(i,j);
            %             newx(i,j)=lu(2,j);
        end
        if newx(i,j)<lu(1,j)
            newx(i,j)=2*lu(1,j)-newx(i,j);
            %             newx(i,j)=lu(1,j);
        end

    end
end


newx=max(newx,repmat(lu(1,1:n),m,1));
newx=min(newx,repmat(lu(2,1:n),m,1));
end