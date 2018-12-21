%% Check order of convergence for 
clear; 
close all; 

%% 

p2 = 2:1:6;                         % 
err2D = ones(size(p2));             % Maxnorm error
tF2 = ones(size(p2));               % Factorization time
tS2 = ones(size(p2));               % Solving time
fill_ratio2 = ones(size(p2));       % fill in ratio 

for p = p2 
    [u, u_ex, err2D(p-1),tF2(p-1),tS2(p-1), fill_ratio2(p-1)] = SolveProblem(p,2,3);
end 

p3 = 2:1:3;                         %
err3D = ones(size(p3));             % Maxnorm error
tF3 = ones(size(p3));               % Factorization time
tS3 = ones(size(p3));               % Solving time
fill_ratio3 = ones(size(p3));       % fill in ratio 

for p= p3
    [u, u_ex, err3D(p-1),tF3(p-1),tS3(p-1),fill_ratio3(p-1)] = SolveProblem(p,3,3);
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
figure; 
plot(p2,tF2)
title('Factorization Time Cholesky Decomposition 2D')
xlabel('p')
ylabel('Time (s)')

figure; 
plot(p3,tF3)
title('Factorization Time Cholesky Decomposition 3D')
xlabel('p')
ylabel('Time (s)')

figure; 
plot(p2,tS2)
title('Factorization Time Solving 2D')
xlabel('p')
ylabel('Time (s)')

figure; 
plot(p3,tS3)
title('Factorization Time Solving 3D')
xlabel('p')
ylabel('Time (s)')

%% Fill ratio analysis
figure; 
plot(p2,tF2)
title('Factorization Time Cholesky Decomposition 2D')
xlabel('p')
ylabel('Time (s)')

figure; 
plot(p3,tF3)
title('Factorization Time Cholesky Decomposition 3D')
xlabel('p')
ylabel('Time (s)')