clear all
close all 

load('14-01-2019/SSOR_mmax200_3D.mat')
load('14-01-2019/SSOR_mmax200_2D.mat')
SSOR_time_2D = tS2;
SSOR_time_3D = tS3;
N_2 = N2;
N_3 = N3;

load('14-01-2019/Cholesky_red_2D.mat')
load('14-01-2019/Cholesky_red_3D.mat')
Cholesky_time_2D = tS2(1:size(N_2,2)) + tF2(1:size(N_2,2));
Cholesky_time_3D = tS3(1:size(N_3,2)) + tF3(1:size(N_3,2));

load('14-01-2019/PCG_mmax200_2D.mat')
load('14-01-2019/PCG_mmax200_3D.mat')
PCG_time_2D = tS2(1:size(N_2,2));
PCG_time_3D = tS3(1:size(N_3,2));

figure;
plot(N_2,Cholesky_time_2D);
hold on;
plot(N_2,SSOR_time_2D);
plot(N_2,PCG_time_2D);
title(['Solving time of all three solvers 2D'])
xlabel('N')
ylabel('time (s)')
set(gca, 'XScale','log')
set(gca, 'YScale', 'log')
legend('Reordered Cholesky', 'SSOR', 'PCG', 'location', 'best')
hold off;

figure;
plot(N_3,Cholesky_time_3D);
hold on;
plot(N_3,SSOR_time_3D);
plot(N_3,PCG_time_3D);
title(['Solving time of all three solvers 3D'])
xlabel('N')
ylabel('time (s)')
set(gca, 'XScale','log')
set(gca, 'YScale', 'log')
legend('Reordered Cholesky', 'SSOR', 'PCG', 'location', 'best')
hold off;