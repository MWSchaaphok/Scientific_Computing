%% Check order of convergence for 
clear; 
close all; 

%% 

p2 = 2:1:2;    
n2 = 2.^p2;                         %
N2 = (n2+ones(size(n2))).^2;
h2 = 1./n2;     
err2D = ones(size(p2));             % Maxnorm error
tF2 = ones(size(p2));               % Factorization time
tS2 = ones(size(p2));               % Solving time
fill_ratio2 = ones(size(p2));       % fill in ratio 
resid2 = ones(size(p2,1),10);      % residual SSOR 
rrf2 = ones(size(p2,1),5);          % residual reduction factor 

for p = p2 
    [u2, u_ex2, err2D(p-1),tF2(p-1),tS2(p-1), fill_ratio2(p-1),resid2(p-1,:),rrf2(p-1,:)] = SolveProblem(p,2,3,'Cholesky',0);
end 


p3 = 2:1:2;                         %
n3 = 2.^p3;                         %
N3 = (n3 + ones(size(n3))).^3;
h3 = 1./n3;
err3D = ones(size(p3));             % Maxnorm error
tF3 = ones(size(p3));               % Factorization time
tS3 = ones(size(p3));               % Solving time
fill_ratio3 = ones(size(p3));       % fill in ratio 
resid3 = ones(size(p3,1),10);        % residual SSOR
rrf3 = ones(size(p3,1),5);          % residual reduction factor 

for p= p3
    [u3, u_ex3, err3D(p-1),tF3(p-1),tS3(p-1),fill_ratio3(p-1), resid3(p-1,:),rrf3(p-1,:)] = SolveProblem(p,3,3,'Cholesky',0);
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
plot(N2,tF2)
title('Factorization Time Cholesky Decomposition 2D')
xlabel('N')
ylabel('Time (s)')

figure; 
plot(N3,tF3)
title('Factorization Time Cholesky Decomposition 3D')
xlabel('N')
ylabel('Time (s)')

figure; 
plot(N2,tS2)
title('Forward/Backward Solving Time 2D')
xlabel('N')
ylabel('Time (s)')

figure; 
plot(N3,tS3)
title('Forward/Backward Solving Time 3D')
xlabel('N')
ylabel('Time (s)')

figure; 
plot(N2,tS2);
title('Solving time SSOR 2D')
xlabel('N')
ylabel('time (s)')
set(gca, 'XScale','log')
set(gca, 'YScale', 'log')
hold on; 
plot(N2, 10^-5*N2.^2); 
legend('t','N^2')
hold off;

%% Fill ratio analysis
figure; 
plot(N2,fill_ratio2)
title('Fill ratio 2D')
xlabel('N')
ylabel('nnz(C)/nnz(A)')

figure; 
plot(N3,fill_ratio3)
title('Fill ratio 3D')
xlabel('N')
ylabel('nnz(C)/nnz(A)')

%% Plot relative residuals SSOR 
figure; 
plot(resid2');
set(gca, 'YScale', 'log')
title('Relative residual ||r||_2/||f||_2 2D')
for i=1:length(n2)
    legendn{i} = sprintf('n=%s',num2str(n2(i)));
end
legend(legendn,'Location','best')
xlabel('m')
ylabel('||r||_2/||f||_2')
%%
% 
figure; 
plot(resid3');
set(gca, 'YScale', 'log')
title('Relative residual ||r||_2/||f||_2 3D')
for i=1:length(n3)
    legendn{i} = sprintf('n=%s',num2str(n3(i)));
end
legend(legendn,'Location','best')
xlabel('m')
ylabel('||r||_2/||f||_2')


%latex_table = latex(vpa(sym(rrf2),3))
