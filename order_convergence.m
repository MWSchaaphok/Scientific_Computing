%% Check order of convergence for 
clear; 
close all; 
%% 
err2D = ones(3,1);

for p2 = 2:1:4
    [u, u_ex, err2D(p2-1)] = SolveProblem(p2,2,3);
end 

err3D = ones(3,1); 
for p3 = 2:1:4
    [u, u_ex, err3D(p3-1)] = SolveProblem(p3,3,3);
end 

figure; 
plot([2,3,4],err2D);
title('Maxnorm |u^h - u^h_{ex}| for 2D')
xlabel('p')
ylabel('maxnorm')

figure; 
plot([2,3,4],err3D);
title('Maxnorm |u^h - u^h_{ex}| for 3D')
xlabel('p')
ylabel('maxnorm')
