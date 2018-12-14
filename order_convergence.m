%% Check order of convergence for 
clear; 
close all; 

%% 

p2 = 2:1:4;
err2D = ones(size(p2));
t2 = ones(size(p2));
for p = p2 
    [u, u_ex, err2D(p-1),t2(p-1)] = SolveProblem(p,2,3);
end 

p3 = 2:1:4; 
err3D = ones(size(p3)); 
t3 = ones(size(p3));
for p= p3
    [u, u_ex, err3D(p-1),t3(p-1)] = SolveProblem(p,3,3);
end 

figure; 
plot(p2,err2D);
title('Maxnorm |u^h - u^h_{ex}| for 2D')
xlabel('p')
ylabel('maxnorm')

figure; 
plot(p3,err3D);
title('Maxnorm |u^h - u^h_{ex}| for 3D')
xlabel('p')
ylabel('maxnorm')

%% Factorization and solving time as function of problem size 
