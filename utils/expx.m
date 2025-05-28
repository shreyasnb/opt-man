function e = expx(x,v)
    a = cos(norm(v))*x + sin(norm(v))*v/norm(v);
    e = a/norm(a);
end