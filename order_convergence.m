%% Check order of convergence for 
clear; 
close all; 

%% 

p2 = 2:1:7;    
n2 = 2.^p2;                         % 
h2 = 1./n2;     
err2D = ones(size(p2));             % Maxnorm error
tF2 = ones(size(p2));               % Factorization time
tS2 = ones(size(p2));               % Solving time
fill_ratio2 = ones(size(p2));       % fill in ratio 

for p = p2 
    [u, u_ex, err2D(p-1),tF2(p-1),tS2(p-1), fill_ratio2(p-1)] = SolveProblem(p,2,3,'Cholesky',0);
end 


p3 = 2:1:2;                         %
n3 = 2.^p3;
h3 = 1./n3;
err3D = ones(size(p3));             % Maxnorm error
tF3 = ones(size(p3));               % Factorization time
tS3 = ones(size(p3));               % Solving time
fill_ratio3 = ones(size(p3));       % fill in ratio 

for p= p3
    [u, u_ex, err3D(p-1),tF3(p-1),tS3(p-1),fill_ratio3(p-1)] = SolveProblem(p,3,3,'Choleksy',0);
end 
 

%%
figure; 
plot(h2,err2D);
title('Maxnorm |u^h - u^h_{ex}| for 2D')
xlabel('h')
ylabel('maxnorm')
set(gca, 'XScale','log')
set(gca, 'YScale', 'log')
hold on; 
plot(h2, h2.^2); 
legend('|e|_2','h^2')
hold off; 

figure; 
set(gca, 'YScale', 'log')
plot(h3,err3D);
title('Maxnorm |u^h - u^h_{ex}| for 3D')
xlabel('h')
ylabel('maxnorm')
set(gca, 'XScale','log')
set(gca, 'YScale', 'log')
hold on; 
plot(h3, h3.^2); 
legend('|e|_2','h^2')
hold off;
%% Factorization and solving time as function of problem size 
figure; 
plot(n2,tF2)
title('Factorization Time Cholesky Decomposition 2D')
xlabel('n')
ylabel('Time (s)')

figure; 
plot(n3,tF3)
title('Factorization Time Cholesky Decomposition 3D')
xlabel('n')
ylabel('Time (s)')

figure; 
plot(n2,tS2)
title('Forward/Backward Solving Time 2D')
xlabel('n')
ylabel('Time (s)')

figure; 
plot(n3,tS3)
title('Forward/Backward Solving Time 3D')
xlabel('n')
ylabel('Time (s)')

%% Fill ratio analysis
figure; 
plot(n2,fill_ratio2)
title('Fill ratio 2D')
xlabel('n')
ylabel('nnz(C)/nnz(A)')

figure; 
plot(n3,fill_ratio3)
title('Fill ratio 3D')
xlabel('n')
ylabel('nnz(C)/nnz(A)')