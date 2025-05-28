function [b, x0, Y0, Y_hat] = init(problem)
    n = problem.dim1;
    k = problem.dim2;
    M = randn(n);
    [Y0,~] = eigs((M+M')/2,k);
    theta = problem.params.theta;
    R = [-theta theta];
    theta0 = rand(1)*(max(R)-min(R)) + min(R);

    if(n==2)
        b = [cos(problem.data.theta_b); sin(problem.data.theta_b)];
        x0 = [cos(theta0); sin(theta0)];
        Y_hat = [1; 0];
    else
        b = [cos(problem.data.theta_b); sin(problem.data.theta_b); 0];
        x0 = [cos(theta0); sin(theta0); 0];
        Y_hat = [1; 0; 0];
    end
end