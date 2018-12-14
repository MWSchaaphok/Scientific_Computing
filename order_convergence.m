%% Check order of convergence
err2D = ones(8,1);

for p = [2,3,4]
    [u, u_ex, err2D(p)] = SolveProblem(p,2,3);
end 

for p = [2,3,4]
    [u, u_ex, err2D(p)] = SolveProblem(p,2,3);
end 
