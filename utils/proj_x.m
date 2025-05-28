function p = proj_x(x,v)
    n = length(x);
    p = (eye(n)-x*x')*v;
end