function problem = init(problem)
    n = problem.dim(1);
    k = problem.dim(2);
    M = randn(n);
    [problem.Y,~] = eigs((M+M')/2,k);
    theta = asin(problem.params.rho);
    R = [-theta theta];
    thetax = rand(1)*(max(R)-min(R)) + min(R);
    theta_0 = problem.data.theta_0;

    if(n==2 && k==1)
        problem.data.b = [cos(problem.data.theta_b); sin(problem.data.theta_b)];
        problem.x0 = [cos(thetax); sin(thetax)];
        problem.data.Y_hat = [cos(theta_0); sin(theta_0)];
    elseif(n==3 && k==1)
        problem.data.b = [cos(problem.data.theta_b); 0; sin(problem.data.theta_b)];
        problem.x0 = [cos(thetax); 0; sin(thetax)];
        problem.data.Y_hat = [1; 0; 0];
    else
        problem.data.b = rand(n,1);
        problem.x0 = [cos(thetax); zeros(n-2,1); sin(thetax)];
        problem.data.Y_hat = eye(n,k);
    end
end