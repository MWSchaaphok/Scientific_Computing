clear all
close all 

load('PCG_RelRes.mat')
PCG_resid2 = resid2;
PCG_resid3 = resid3;
load('SSOR_RelRes.mat')
SSOR_resid2 = resid2;
SSOR_resid3 = resid3;

figure;
plot(PCG_resid2');
hold on
plot(SSOR_resid2');
set(gca, 'YScale', 'log')
title(['Relative residuals ||r||_2/||f||_2 for SSOR and PCG 2D'])
for i=1:length(n2)
    legendn{i} = sprintf('PCG, n=%s',num2str(n2(i)));
end
for i=1:length(n2)
    legendn{i+length(n2)} = sprintf('SSOR, n=%s',num2str(n2(i)));
end
legend(legendn,'Location','best')
xlabel('m')
ylabel('||r||_2/||f||_2')


%% 3D %%
figure;
plot(PCG_resid3');
hold on
plot(SSOR_resid3');
set(gca, 'YScale', 'log')
title(['Relative residuals ||r||_2/||f||_2 for SSOR and PCG 2D'])
for i=1:length(n2)
    legendn{i} = sprintf('PCG, n=%s',num2str(n2(i)));
end
for i=1:length(n2)
    legendn{i+length(n2)} = sprintf('SSOR, n=%s',num2str(n2(i)));
end
legend(legendn,'Location','best')
xlabel('m')
ylabel('||r||_2/||f||_2')

