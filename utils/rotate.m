function [xr,yr] = rotate(x, y, theta)
    n = length(x);
    xmin = min(x);
    xmax = max(x);
    xr = x.*cos(theta) - y.*sin(theta);
    yr = x.*sin(theta) + y.*cos(theta);
end