function y=repair(x,bounds)
    lower=repmat(bounds(1,:),size(x,1),1);
    upper=repmat(bounds(2,:),size(x,1),1);
    
    y=x;
    
    % reflected from the bounds
    pos= x<lower;
    y(pos)=2*lower(pos) -x(pos);  % symetric at x0=l
    
    pos= x>upper;
    y(pos)=2*upper(pos) -x(pos);  % symetric at x0=l
    
    % make sure solutions are in the bounds
    y=max(y, lower);
    y=min(y, upper);
end